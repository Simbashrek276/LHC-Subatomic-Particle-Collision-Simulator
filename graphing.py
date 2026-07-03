# Data analysis stage: read the events that simulation.py logged and draw an
# energy histogram for each of the four particles in the 2-photon/2-proton
# events. Run simulation.py first so that no_higgs_events.txt exists.

import matplotlib.pyplot as plt

DATA_FILE = "no_higgs_events.txt"


def read_energies(filename):
    """Pull every particle energy out of the data file, in the order written.

    Each logged event is four particle rows, so the energies come back grouped
    as [photon1, photon2, proton1, proton2] for one event, then the next, and
    so on -- exactly the order the friend's drawing code below expects.
    """
    all_energies = []
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
            ):
                continue
            # A particle row looks like "photon   12.345678   ...", so the
            # second column is the energy.
            all_energies.append(float(line.split()[1]))
    return all_energies


all_energies = read_energies(DATA_FILE)

# --- Original drawing code (kept in the same structure as before) ----------
# Graph all energies from the logged events. The energies arrive as one long
# list, so we deal them out into four buckets: every 4th value belongs to the
# same particle position (photon 1, photon 2, proton 1, proton 2).
particle_energies = [[], [], [], []]

for i, energy in enumerate(all_energies):
    particle_energies[i % 4].append(energy)

particle_names = [
    "Photon 1",
    "Photon 2",
    "Proton 1",
    "Proton 2"
]

# Draw and save one energy histogram per particle position.
for i in range(4):
    plt.figure(figsize=(8, 5))
    plt.hist(particle_energies[i], bins=20)

    plt.xlabel(f"Energy of {particle_names[i]} Particle (TeV)")
    plt.ylabel("Number of Events")
    plt.title(f"Energy Distribution of {particle_names[i]} Particle")

    plt.tight_layout()
    filename = (
        particle_names[i]
        .lower()
        .replace(" ", "_")
        + "_energy_histogram.png"
    )
    plt.savefig(filename)
    print(f"Histogram saved to {filename}")
