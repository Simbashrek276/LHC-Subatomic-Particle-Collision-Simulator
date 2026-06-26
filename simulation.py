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

        n = len(particles)
        energies = []
        angles = []
        px_list = []
        pz_list = []
        total_px = 0.0
        total_pz = 0.0

        REMAINING_ENERGY = TOTAL_ENERGY

        # Kinematics for the first n-1 particles
        for j in range(n - 1):
            E = random.uniform(0.0, REMAINING_ENERGY)
            REMAINING_ENERGY -= E
            theta = random.uniform(0, 2 * PI)

            pz = E * math.cos(theta)
            px = E * math.sin(theta)

            energies.append(E)
            angles.append(theta)
            px_list.append(px)
            pz_list.append(pz)

            total_px += px
            total_pz += pz

        # Kinematics for the last particle
        E_last = REMAINING_ENERGY
        px_last = -total_px
        pz_last = -total_pz
        theta_last = math.atan2(px_last, pz_last)

        energies.append(E_last)
        angles.append(theta_last)
        px_list.append(px_last)
        pz_list.append(pz_last)

        # Does the whole event pass the energy + angle cuts?
        passes_cuts = all(
            E >= MIN_ENERGY and is_valid_angle(theta)
            for E, theta in zip(energies, angles)
        )

        # Build the particle rows once (shared by console and file).
        body = [
            f"State: {state}\n",
            f"N_Particles: {n}\n",
            f"{'Particle':<12} {'Energy(TeV)':<15} {'Angle(rad)':<15} {'p_x(TeV)':<15} {'p_z(TeV)':<15}\n",
        ]
        for particle, E, theta, px, pz in zip(
            particles, energies, angles, px_list, pz_list
        ):
            body.append(
                f"{particle:<12} {E:<15.6f} {theta % (2*PI):<15.6f} {px:<15.6f} {pz:<15.6f}\n"
            )

        # Print EVERY event/combination to the console, tagged with whether it
        # passed the kinematic cuts.
        cut_tag = "PASSED CUTS" if passes_cuts else "FAILED CUTS"
        print(f"___COLLISION: {total_simulated_collisions}___ [{cut_tag}]")
        for line in body:
            print(line, end="")
        print()

        # Export ONLY the 2-photon/2-proton states (Higgs + No Higgs) that pass
        # the cuts, up to the target number of logged events.
        if state in ("Higgs", "No Higgs") and passes_cuts:
            logged_count += 1
            f.write(f"___EVENT_ID: {logged_count}___\n")
            for line in body:
                f.write(line)