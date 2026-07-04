import matplotlib.pyplot as plt
import math

DATA_FILE = "no_higgs_events.txt"


def read_angles(filename):
    """Pull every particle energy out of the data file, in the order written.

    Each logged event is four particle rows, so the energies come back grouped
    as [photon1, photon2, proton1, proton2] for one event, then the next, and
    so on -- exactly the order the friend's drawing code below expects.
    """
    all_angles = []
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
            ):
                continue
            # A particle row looks like "photon   12.345678   ...", so the
            # second column is the energy.
            all_angles.append(float(line.split()[2]))
    return all_angles


all_angles = read_angles(DATA_FILE)

# --- Original drawing code (kept in the same structure as before) ----------
# Graph all angles from the logged events. The angles arrive as one long
# list, so we deal them out into four buckets: every 4th value belongs to the
# same particle position (photon 1, photon 2, proton 1, proton 2).
particle_angles = [[], [], [], []]

for i, angle in enumerate(all_angles):
    particle_angles[i % 4].append(math.cos(angle))

particle_names = [
    "Photon 1",
    "Photon 2",
    "Proton 1",
    "Proton 2"
]

# Draw and save one energy histogram per particle position.
for i in range(4):
    plt.figure(figsize=(8, 5))
    plt.hist(particle_angles[i], bins=20)

    plt.xlabel(f"Cos of Angle of {particle_names[i]} Particle (rad)")
    plt.ylabel("Number of Events")
    plt.title(f"Cos of Angle Distribution of {particle_names[i]} Particle")

    plt.tight_layout()
    filename = (
        particle_names[i]
        .lower()
        .replace(" ", "_")
        + "_cos_angle_histogram.png"
    )
    plt.savefig(filename)
    print(f"Histogram saved to {filename}")