import matplotlib.pyplot as plt

#draw the graph such that instead of coloring the histogram bars, 
#we only draw the top line of the bars. We stack 3 graphs on top of each 
#other for 1 graph. 
#under the graph should be a ratio graph that shows the ratio
#between n4/n3, n5/n3, and n5/n4.

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


energies = read_2to3(DATA_FILE, column=1)
total_events = len(energies[0])   # one entry per 2 to 3 event
print("Total 2 to 3 events:", total_events)


#historgram for particle E45 only but by using 13.6 - E3
plt.figure(figsize=(8, 5))
E45 = []
for i in range(len(energies[0])):
    E45.append(13.6 - energies[0][i])   # 13.6 TeV is the full energy
plt.hist(E45, bins=40, range=(0, 13.6), alpha=0.5, label=PARTICLE_LABELS[1])   # 13.6 TeV is the full energy

plt.xlabel(f"Energy of {PARTICLE_LABELS[1]} (TeV)")
plt.ylabel("Number of events")
plt.title(f"Energy Distribution for the 2 to 3 Case (total events = {total_events})")

plt.tight_layout()
filename = f"2to3_energy_particle (45).png"   # particle 4 and 5
plt.savefig(filename)
print("Histogram saved to", filename)
