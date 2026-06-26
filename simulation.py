import random
import math

TOTAL_ENERGY = 13.6
N_COLLISIONS = 10

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
        return ["photon", "photon", "proton", "proton"], None
    elif r < 0.50:
        return ["positron", "electron", "proton", "proton"], None
    elif r < 0.68:
        return ["muon", "antimuon", "proton", "proton"], None
    else:
        return ["neutron", "antineutron", "proton", "proton"], None

for i in range(N_COLLISIONS):
    particles, state = choose_final_state()
    n = len(particles)
    energies = []
    angles = []
    px_list = []
    pz_list = []
    total_px = 0.0
    total_pz = 0.0

    REMAINING_ENERGY = TOTAL_ENERGY
    for j in range(n - 1):
        E = random.uniform(0.0, REMAINING_ENERGY) # Fix this since 
        #conservation of energy is not conserved. After every
        #random, you have to minus the first randomed energy
        #then continue random the remaining particles' energy
        REMAINING_ENERGY -= E
        theta = random.uniform(0, 2 * math.pi)

        #z la truc ngang, x la truc doc
        pz = E * math.cos(theta)
        px = E * math.sin(theta)

        energies.append(E)
        angles.append(theta)
        px_list.append(px)
        pz_list.append(pz)

        total_px += px
        total_pz += pz

    px_last = -total_px
    pz_last = -total_pz
    # Last particle gets the leftover energy so sum(energies) == TOTAL_ENERGY.
    E_last = REMAINING_ENERGY
    energies.append(E_last)
    angles.append(math.atan2(px_last, pz_last))

    px_list.append(px_last)
    pz_list.append(pz_last)

    print(f"\nCOLLISION {i + 1}")
    if state is not None:
        print(f"State: {state}")
    print("Final state is:", " + ".join(particles))
    for particle, E, theta, px, pz in zip(
        particles, energies, angles, px_list, pz_list
    ):
        print(f"Particle: {particle}")
        print(f"___Energy (TeV): {E}")
        print(f"___Angle (rad): {theta}")
        print(f"___p_x (TeV): {px}")
        print(f"___p_z (TeV): {pz}")