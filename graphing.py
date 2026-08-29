import matplotlib.pyplot as plt

DATA_FILE = "events.txt"


def read_2to2_angles(filename):
    # collect 2to2 collisions' angles
    angles = []
    n_particles = 0
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith("N_Particles"):
                n_particles = int(line.split(":")[1])
            elif line and not line.startswith("___") and not line.startswith("Particle"):
                # this is a particle row, keep it only if the event has 2 particles
                if n_particles == 2:
                    parts = line.split()
                    angles.append(float(parts[2]))   # the 3rd column is the angle in degrees
    return angles


angles = read_2to2_angles(DATA_FILE)
print("Number of 2 to 2 particles:", len(angles))

plt.figure(figsize=(8, 5))
plt.hist(angles, bins=36, range=(0, 360))   # 36 bins means 10 degrees each

plt.xlabel("Angle (degrees)")
plt.ylabel("Number of Particles")
plt.title("Angle Distribution for the 2 to 2 Case")

plt.tight_layout()
plt.savefig("2to2_angle_histogram.png")
print("Histogram saved to 2to2_angle_histogram.png")
