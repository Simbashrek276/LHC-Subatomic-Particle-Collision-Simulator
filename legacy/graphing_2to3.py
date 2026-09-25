import matplotlib.pyplot as plt

DATA_FILE = "events.txt"

# In the 2 to 3 case the three particles are a photon (particle 3) and two
# protons (particle 4 and particle 5). We keep them apart so each gets its own
# histogram.
PARTICLE_LABELS = ["Photon (particle 3)", "Proton (particle 4)", "Proton (particle 5)"]


def read_2to3(filename, column):
    # Return one list per particle for the 3 particle events.
    # column 1 is the energy and column 2 is the angle.
    per_particle = [[], [], []]
    n_particles = 0
    index = 0
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith("N_Particles"):
                n_particles = int(line.split(":")[1])
                index = 0
            elif line and not line.startswith("___") and not line.startswith("Particle"):
                # this is a particle row, keep it only if the event has 3 particles
                if n_particles == 3:
                    parts = line.split()
                    per_particle[index].append(float(parts[column]))
                    index += 1
    return per_particle


angles = read_2to3(DATA_FILE, column=2)
total_events = len(angles[0])   # one entry per 2 to 3 event
print("Total 2 to 3 events:", total_events)

# One angle histogram per particle.
for i in range(3):
    plt.figure(figsize=(8, 5))
    plt.hist(angles[i], bins=36, range=(0, 360))   # 36 bins means 10 degrees each

    plt.xlabel(f"Flight angle of {PARTICLE_LABELS[i]} (degrees)")
    plt.ylabel("Number of events")
    plt.title(f"Angle Distribution for the 2 to 3 Case (total events = {total_events})")

    plt.tight_layout()
    filename = f"2to3_angle_particle{i + 3}.png"   # particle 3, 4 and 5
    plt.savefig(filename)
    print("Histogram saved to", filename)
