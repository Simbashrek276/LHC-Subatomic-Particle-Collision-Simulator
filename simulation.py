# I edited the file to 2 functions: 
#
#   1. collision()  -- throws the dice for a single collision: picks a final
#                       state and shares out energy and momentum among the
#                       particles.
#   2. detector()   -- plays the role of the detector: it applies the energy
#                       and angle cuts, prints what it saw, and records the
#                       events we care about into no_higgs_events.txt.
#
# Run this file first to produce the data file, then run graphing.py to draw the histograms thing.

import math
import random

TOTAL_ENERGY = 13.6  # TeV
BEAM_ENERGY = 6.8  # TeV -- each of the two incoming protons carries half the total
TARGET_LOGGED_EVENTS = 4000
OUTPUT_FILE = "no_higgs_events.txt"

MIN_ENERGY = 0.02 
PI = math.pi
ANGLE_RANGES = [
    (1 / 18 * PI, 17 / 18 * PI),   # 10 to 170 degrees
    (19 / 18 * PI, 35 / 18 * PI),  # 190 to 350 degrees
]


# Each final state returns (particles, state_label). "Higgs" is only a STATE
# label for the photon/photon/proton/proton final state
def choose_final_state():
    r = random.random()
    if r < 0.1: #past is 0.25
        return ["photon", "proton", "proton"], None
    elif r < 0.4: #past is 0.26
        return ["photon", "photon", "proton", "proton"], "Higgs"
    elif r < 0.8: #past is 0.30
        return ["photon", "photon", "proton", "proton"], "No Higgs"
    elif r < 0.9: #past is 0.50
        return ["positron", "electron", "proton", "proton"], None
    elif r < 0.95: #past is 0.68
        return ["muon", "antimuon", "proton", "proton"], None
    else:
        return ["neutron", "antineutron", "proton", "proton"], None


def is_valid_angle(theta):
    t = theta % (2 * PI)  # normalize to [0, 2*pi) first
    for low, high in ANGLE_RANGES:
        if low <= t <= high:
            return True
    return False


def beam_split_energies():
    """Work out four energies by asking how much of itself each proton gave up.

    Picture the actual smash: proton 1 and proton 2 come in carrying 6.8 TeV
    each. Each one hands a slice of its energy over to a photon and keeps the
    rest. So we only need to roll two numbers -- x1 for the first proton and x2
    for the second -- and everything else follows:

        photon 3 = x1 * 6.8        proton 5 = (1 - x1) * 6.8
        photon 4 = x2 * 6.8        proton 6 = (1 - x2) * 6.8

    What a proton gives away plus what it keeps is just its whole 6.8 TeV, so
    the four energies always add straight back up to the full 13.6 TeV. 
    """
    x1 = random.random()
    x2 = random.random()
    return [
        x1 * BEAM_ENERGY,        # photon 3 -- the slice proton 1 gave up
        x2 * BEAM_ENERGY,        # photon 4 -- the slice proton 2 gave up
        (1 - x1) * BEAM_ENERGY,  # proton 5 -- what proton 1 held on to
        (1 - x2) * BEAM_ENERGY,  # proton 6 -- what proton 2 held on to
    ]


def check_energy_adds_up(energies):
    total = sum(energies)
    if not math.isclose(total, TOTAL_ENERGY, abs_tol=1e-9):
        raise ValueError(
            f"Energy does not add up: particles total {total:.9f} TeV, "
            f"but the collision started with {TOTAL_ENERGY} TeV"
        )
    return total


def collision():
    """Simulate a single collision and return everything it produced.

    We pick a final state, then hand out energy and momentum. There are two
    ways the energy gets shared, depending on what came out:

      * Two-photon states use beam_split_energies(): each incoming proton gives
        a slice of its 6.8 TeV to a photon and keeps the rest.
      * Everything else falls back to the older approach, where each particle
        takes a random bite out of whatever budget is still on the table.

    Either way the energies add up to TOTAL_ENERGY. The last particle's
    momentum is then set to balance the others so the momenta sum to zero.
    """
    particles, state = choose_final_state()
    n = len(particles)

    energies = []
    angles = []
    px_list = []
    pz_list = []
    total_px = 0.0
    total_pz = 0.0

    remaining_energy = TOTAL_ENERGY

    # The two-photon states get their energies from the beam-splitting picture
    # above. Every other final state still shares out the budget the old way.
    planned_energies = None
    if state in ("Higgs", "No Higgs"):
        planned_energies = beam_split_energies()

    # Share energy/momentum among the first n-1 particles.
    for j in range(n - 1):
        if planned_energies is not None:
            E = planned_energies[j]
        else:
            E = random.uniform(0.0, remaining_energy)
        remaining_energy -= E
        theta = random.uniform(0, 2 * PI)

        pz = E * math.cos(theta)
        px = E * math.sin(theta)

        energies.append(E)
        angles.append(theta)
        px_list.append(px)
        pz_list.append(pz)

        total_px += px
        total_pz += pz

    # The last particle gets momentum that cancels everyone else's. Its energy is
    # proton 6's share in the two-photon case, or the leftover budget otherwise.
    if planned_energies is not None:
        E_last = planned_energies[n - 1]
    else:
        E_last = remaining_energy
    px_last = -total_px
    pz_last = -total_pz
    theta_last = math.atan2(px_last, pz_last)

    energies.append(E_last)
    angles.append(theta_last)
    px_list.append(px_last)
    pz_list.append(pz_last)

    # Now that everyone has their share, make sure it still adds back up to 13.6.
    check_energy_adds_up(energies)

    return {
        "particles": particles,
        "state": state,
        "energies": energies,
        "angles": angles,
        "px_list": px_list,
        "pz_list": pz_list,
    }


def detector(event, f, collision_number, logged_so_far):
    """Ham nay se look at one collision, apply the cuts, and record it if we want it.

    Every collision is printed to the console (tagged with whether it passed
    the cuts). Only the 2-photon/2-proton states (Higgs and No Higgs) that pass
    all cuts get written to the output file. Returns the event's energies when
    it was logged, otherwise None.
    """
    particles = event["particles"]
    state = event["state"]
    energies = event["energies"]
    angles = event["angles"]
    px_list = event["px_list"]
    pz_list = event["pz_list"]
    n = len(particles)

    # An event passes only if EVERY particle clears the energy and angle cuts.
    passes_cuts = all(
        E >= MIN_ENERGY and is_valid_angle(theta)
        for E, theta in zip(energies, angles)
    )

    # Build the particle rows once; both the console and the file reuse them.
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

    # Show every collision we simulated, passed or failed. The running total is
    # printed to screen only -- the data file keeps its original layout so the
    # graphing scripts can still read it.
    total_energy = sum(energies)
    cut_tag = "PASSED CUTS" if passes_cuts else "FAILED CUTS"
    print(
        f"___COLLISION: {collision_number}___ [{cut_tag}] "
        f"[Total Energy: {total_energy:.6f} TeV]"
    )
    for line in body:
        print(line, end="")
    print()

    # i keep only the 2-photon/2-proton events (Higgs + No Higgs) that pass cuts.
    if state in ("Higgs", "No Higgs") and passes_cuts:
        event_id = logged_so_far + 1
        f.write(f"___EVENT_ID: {event_id}___\n")
        for line in body:
            f.write(line)
        return energies

    return None


def main():
    logged_count = 0
    total_collisions = 0

    # I let the simulation run 
    # until there are enough number of wanted 2photon2proton events.
    with open(OUTPUT_FILE, "w") as f:
        while logged_count < TARGET_LOGGED_EVENTS:
            total_collisions += 1
            event = collision()
            energies = detector(event, f, total_collisions, logged_count)
            if energies is not None:
                logged_count += 1


if __name__ == "__main__":
    main()
