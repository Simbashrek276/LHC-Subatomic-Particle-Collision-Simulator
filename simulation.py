import math
import random

TOTAL_ENERGY = 13.6
TARGET_LOGGED_EVENTS = 100
OUTPUT_FILE = "events.txt"
MIN_ENERGY = 0.02

PI = math.pi

ANGLE_RANGES = [
    (10 * PI / 180, 170 * PI / 180),
    (190 * PI / 180, 350 * PI / 180),
]

MASS = {
    "photon": 0.0,
    "proton": 0.938,
    "positron": 0.000511,
    "electron": 0.000511,
    "muon": 0.1057,
    "antimuon": 0.1057,
    "neutron": 0.9396,
    "antineutron": 0.9396
}

def choose_final_state():
    r = random.random()
    if r < 0.1:
        return ["photon", "proton", "proton"]
    elif r < 0.4:
        return ["photon", "photon", "proton", "proton"]
    elif r < 0.8:
        return ["photon", "photon", "proton", "proton"]
    elif r < 0.9:
        return ["positron", "electron", "proton", "proton"]
    elif r < 0.95:
        return ["muon", "antimuon", "proton", "proton"]
    else:
        return ["neutron", "antineutron", "proton", "proton"]

def is_valid_angle(theta):
    theta = theta % (2 * PI)
    for low, high in ANGLE_RANGES:
        if low <= theta <= high:
            return True
    return False

def solve_scale(vectors, masses):
    # This function finds a scale factor that makes the total
    # energy of all particles equal to TOTAL_ENERGY (13.6 TeV).
    #
    # The original momentum vectors are randomly generated and
    # do not necessarily have enough energy. Instead of changing
    # their directions, we multiply all momenta by the same
    # scale factor.

    def total_energy(scale):
        # Calculate the total energy of all particles
        # after multiplying their momenta by the scale factor.
        total = 0.0

        for vector, mass in zip(vectors, masses):

            # Scale the x and z components of momentum.
            px = scale * vector[0]
            pz = scale * vector[1]

            # Calculate the magnitude of the momentum:
            # p = sqrt(px² + pz²)
            p = math.sqrt(px * px + pz * pz)

            # Calculate the relativistic energy:
            # E = sqrt(p² + m²)
            #
            # The result is added to the total energy
            # of all particles.
            total += math.sqrt(p * p + mass * mass)

        return total

    # Start by searching for a scale factor between 0 and 1.
    low = 0.0
    high = 1.0

    # If scale = 1 does not provide enough total energy,
    # keep doubling the upper limit until the total energy
    # is greater than or equal to 13.6 TeV.
    #
    # Example:
    # scale = 1  → energy too low
    # scale = 2  → energy too low
    # scale = 4  → energy high enough
    #
    # Now the correct scale must be somewhere between 2 and 4.
    while total_energy(high) < TOTAL_ENERGY:
        high *= 2.0

    # Use binary search to find the scale factor that gives
    # a total energy as close as possible to 13.6 TeV.
    #
    # Each iteration cuts the possible range in half.
    for _ in range(100):

        # Try the value halfway between low and high.
        middle = (low + high) / 2.0

        # If this scale produces too little energy,
        # the correct scale must be larger.
        if total_energy(middle) < TOTAL_ENERGY:
            low = middle

        # Otherwise, the scale is large enough, so the
        # correct value must be at or below this value.
        else:
            high = middle
    # low and high are now extremely close to the correct
    # scale factor, so return their midpoint.
    return (low + high) / 2.0

def generate_final_state(particles):
    masses = [MASS[p] for p in particles]
    n = len(particles)
    if sum(masses) > TOTAL_ENERGY:
        return None
    vectors = []
    for _ in range(n - 1): #random momenta
        theta = random.uniform(0, 2 * PI)
        magnitude = random.uniform(0.1, 1.0)
        px = magnitude * math.sin(theta)
        pz = magnitude * math.cos(theta)
        vectors.append((px, pz))
    total_px = sum(v[0] for v in vectors)
    total_pz = sum(v[1] for v in vectors)
    vectors.append((-total_px, -total_pz))
    scale = solve_scale(vectors, masses)

    energies = []
    px_list = []
    pz_list = []
    angles = []

    for vector, mass in zip(vectors, masses):
        px = scale * vector[0]
        pz = scale * vector[1]

        momentum = math.sqrt(px * px + pz * pz)
        energy = math.sqrt(momentum * momentum + mass * mass)

        theta = math.atan2(px, pz)

        if theta < 0:
            theta += 2 * PI

        px_list.append(px)
        pz_list.append(pz)
        energies.append(energy)
        angles.append(theta)

    return {
        "particles": particles,
        "energies": energies,
        "angles": angles,
        "px_list": px_list,
        "pz_list": pz_list
    }

def two_to_three(particles):
    if len(particles) != 3:
        return None

    return generate_final_state(particles)

def two_to_four(particles):
    if len(particles) != 4:
        return None

    return generate_final_state(particles)

def collision():
    particles = choose_final_state()
    if len(particles) == 3:
        return two_to_three(particles)
    elif len(particles) == 4:
        return two_to_four(particles)
    return None

def detector(event, f, collision_number, logged_so_far):
    if event is None:
        return None
    particles = event["particles"]
    energies = event["energies"]
    angles = event["angles"]
    px_list = event["px_list"]
    pz_list = event["pz_list"]
    passes_cuts = all(
        E >= MIN_ENERGY and is_valid_angle(theta)
        for E, theta in zip(energies, angles)
    )
    total_energy = sum(energies)
    total_px = sum(px_list)
    total_pz = sum(pz_list)
    energy_ok = math.isclose(
        total_energy,
        TOTAL_ENERGY,
        abs_tol=1e-8
    )
    momentum_ok = (
        math.isclose(total_px, 0.0, abs_tol=1e-8)
        and
        math.isclose(total_pz, 0.0, abs_tol=1e-8)
    )
    if not energy_ok or not momentum_ok:
        return None
    if not passes_cuts:
        return None
    event_id = logged_so_far + 1
    f.write(
        f"___EVENT_ID: {event_id}___\n"
    )
    f.write(
        f"N_Particles: {len(particles)}\n"
    )
    f.write(
        f"{'Particle':<15}"
        f"{'Energy(TeV)':<15}"
        f"{'Angle(deg)':<15}"
        f"{'p_x(TeV)':<15}"
        f"{'p_z(TeV)':<15}\n"
    )
    for particle, E, theta, px, pz in zip(
        particles,
        energies,
        angles,
        px_list,
        pz_list
    ):
        angle_degrees = theta * 180 / PI
        f.write(
            f"{particle:<15}"
            f"{E:<15.6f}"
            f"{angle_degrees:<15.6f}"
            f"{px:<15.6f}"
            f"{pz:<15.6f}\n"
        )
    f.write("\n")
    return energies

logged_count = 0
total_collisions = 0

with open(OUTPUT_FILE, "w") as f:
    for collision_number in range(10000000):
        if logged_count >= TARGET_LOGGED_EVENTS:
            break
        total_collisions += 1
        event = collision()
        energies = detector(
            event,
            f,
            total_collisions,
            logged_count
        )
        if energies is not None:
            logged_count += 1

print("Simulation finished.")
print(f"Total collisions: {total_collisions}")
print(f"Accepted events: {logged_count}")
print(f"Saved to: {OUTPUT_FILE}")