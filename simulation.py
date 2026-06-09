import random
import math

TOTAL_ENERGY = 13.6
N_COLLISIONS = 10000

def choose_final_state():
    r = random.random()
    if r < 0.25:
        return ["photon", "proton", "proton"]
    elif r < 0.26:
        return ["photon", "photon", "proton", "proton", "Higgs"]
    elif r < 0.30:
        return ["photon", "photon", "proton", "proton"]
    elif r < 0.50:
        return ["positron", "electron", "proton", "proton"]
    elif r < 0.68:
        return ["muon", "antimuon", "proton", "proton"]
    else:
        return ["neutron", "antineutron", "proton", "proton"]
for i in range(N_COLLISIONS):
    particles = choose_final_state()
    n = len(particles)
    energies = []
    angles = []
    px_list = []
    pz_list = []
    total_px = 0.0
    total_pz = 0.0
    for j in range(n - 1):
        E = random.uniform(0.0, TOTAL_ENERGY)
        theta = random.uniform(0, 2 * math.pi)
        px = E * math.cos(theta)
        pz = E * math.sin(theta)
        energies.append(E)
        angles.append(theta)
        px_list.append(px)
        pz_list.append(pz)
        total_px += px
        total_pz += pz
    px_last = -total_px
    pz_last = -total_pz
    E_last = math.sqrt(px_last**2 + pz_last**2)
    energies.append(E_last)
    angles.append(math.atan2(pz_last, px_last))
    px_list.append(px_last)
    pz_list.append(pz_last)
    print(f"\nCollision {i + 1}")
    print("Final state is:", " + ".join(particles))
    for particle, E, theta, px, pz in zip(
        particles, energies, angles, px_list, pz_list
    ):
        print(f"Particle: {particle}")
        print(f"Energy (TeV): {E}")
        print(f"Angle (rad): {theta}")
        print(f"p_x (TeV): {px}")
        print(f"p_z (TeV): {pz}")