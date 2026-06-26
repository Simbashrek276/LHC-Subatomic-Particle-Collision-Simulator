import math
import random

TOTAL_ENERGY = 13.6  # TeV
TARGET_LOGGED_EVENTS = 5
OUTPUT_FILE = "no_higgs_events.txt"

# Cuts
MIN_ENERGY = 0.02  # 20 GeV converted to TeV
PI = math.pi
ANGLE_RANGES = [
    (1 / 18 * PI, 17 / 18 * PI),  # 10 to 170 degrees
    (19 / 18 * PI, 35 / 18 * PI),  # 190 to 350 degrees
]


# Each final state returns (particles, state_label). "Higgs" is only a STATE
# label for the photon/photon/proton/proton final state -- it is not a particle,
# so it carries no energy or momentum of its own
def choose_final_state():
    r = random.random()
    if r < 0.25:
        return ["photon", "proton", "proton"], None
    elif r < 0.26:
        return ["photon", "photon", "proton", "proton"], "Higgs"
    elif r < 0.30:
        return ["photon", "photon", "proton", "proton"], "No Higgs"
    elif r < 0.50:
        return ["positron", "electron", "proton", "proton"], None
    elif r < 0.68:
        return ["muon", "antimuon", "proton", "proton"], None
    else:
        return ["neutron", "antineutron", "proton", "proton"], None


def is_valid_angle(theta):
    # Normalize to [0, 2*pi) for validation against the cut windows
    t = theta % (2 * PI)
    for low, high in ANGLE_RANGES:
        if low <= t <= high:
            return True
    return False


logged_count = 0
total_simulated_collisions = 0

with open(OUTPUT_FILE, "w") as f:
    while logged_count < TARGET_LOGGED_EVENTS:
        total_simulated_collisions += 1
        particles, state = choose_final_state()

        # Phase 1 Filter: We only care about logging the "No Higgs" 4-body state
        if state != "No Higgs":
            continue

        n = len(particles)
        energies = []
        angles = []
        px_list = []
        pz_list = []
        total_px = 0.0
        total_pz = 0.0

        REMAINING_ENERGY = TOTAL_ENERGY
        event_is_valid = True

        # Kinematics for the first n-1 particles
        for j in range(n - 1):
            E = random.uniform(0.0, REMAINING_ENERGY)
            REMAINING_ENERGY -= E
            theta = random.uniform(0, 2 * PI)

            # Strict cut verification per particle
            if E < MIN_ENERGY or not is_valid_angle(theta):
                event_is_valid = False
                break

            pz = E * math.cos(theta)
            px = E * math.sin(theta)

            energies.append(E)
            angles.append(theta)
            px_list.append(px)
            pz_list.append(pz)

            total_px += px
            total_pz += pz

        if not event_is_valid:
            continue  # Discard if internal particle fails cuts

        # Kinematics for the last particle
        E_last = REMAINING_ENERGY
        px_last = -total_px
        pz_last = -total_pz
        theta_last = math.atan2(px_last, pz_last)

        # Phase 2 Filter: Check cuts for the final balancing particle
        if E_last < MIN_ENERGY or not is_valid_angle(theta_last):
            continue

        energies.append(E_last)
        angles.append(theta_last)
        px_list.append(px_last)
        pz_list.append(pz_last)

        # If it passes all criteria, officially log it
        logged_count += 1

        f.write(f"EVENT_ID: {logged_count}\n")
        f.write(f"State: {state}\n")
        f.write(f"N_Particles: {n}\n")
        f.write(
            f"{'Particle':<12} {'Energy(TeV)':<15} {'Angle(rad)':<15} {'p_x(TeV)':<15} {'p_z(TeV)':<15}\n"
        )

        for particle, E, theta, px, pz in zip(
            particles, energies, angles, px_list, pz_list
        ):
            f.write(
                f"{particle:<12} {E:<15.6f} {theta % (2*PI):<15.6f} {px:<15.6f} {pz:<15.6f}\n"
            )