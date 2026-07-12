import matplotlib.pyplot as plt
import math

DATA_FILE = "no_higgs_events.txt"


def read_data(filename):
    """Pull every particle energy, px, pz out of the data file, in the order written.

    Each logged event is four particle rows, so the energies come back grouped
    as [photon1, photon2, proton1, proton2] for one event, then the next, and
    so on -- exactly the order the friend's drawing code below expects.
    """
    all_energies = []
    all_px = []
    all_pz = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            # The file also has header/label lines; only the particle rows hold
            # real data, so skip everything else.
            if (
                not line
                or line.startswith("___")
                or line.startswith("State:")
                or line.startswith("N_Particles:")
                or line.startswith("Particle")
                or line.startswith("Energy(TeV)")
                or line.startswith("Angle(rad)")
                or line.startswith("p_x(TeV)")
            ):
                continue

            parts = line.split()
            all_energies.append(float(parts[1]))
            all_px.append(float(parts[3]))
            all_pz.append(float(parts[4]))

    return all_energies, all_px, all_pz


all_energies, all_px, all_pz = read_data(DATA_FILE)

mass_1_2 = []

for i in range(0, len(all_energies), 4):

    E = all_energies[i] + all_energies[i + 1]
    px = all_px[i] + all_px[i + 1]
    pz = all_pz[i] + all_pz[i + 1]

    m2 = E**2 - px**2 - pz**2

    if m2 > 0 and m2<1:
        mass_1_2.append(math.sqrt(m2))

print("Number of events: ", len(mass_1_2))

# draw & save one mass histogram
plt.figure(figsize=(8, 5))
plt.hist(mass_1_2, bins=50)

plt.xlabel(f"Diphoton Invariant Mass (TeV/c^2)")
plt.ylabel("Number of Events")
plt.title(f"Diphoton Invariant Mass Distribution")

plt.tight_layout()
plt.savefig("diphoton_mass_histogram.png")