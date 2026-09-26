"""Shared reader for events.txt. Used by the step-outline graphing scripts.

The old graphing scripts each re-implemented their own line parser and each one
hard coded a column number. That broke when the simulation went to 3D, because
the row layout gained two columns. This module owns the layout in one place so a
future format change is a one line fix.

The 3D row layout written by simulation.write_event is

    token  0        1            2           3          4         5         6
           name     Energy(TeV)  Theta(deg)  Phi(deg)   p_x(TeV)  p_y(TeV)  p_z(TeV)

theta is the polar angle from the beam (z) axis and runs 0 to 180 degrees.
phi is the azimuthal angle around the beam, in the x-y plane, 0 to 360 degrees.

Events are grouped by how many particles came out, because that is what picks
the generator in kinematics.py. Channel 2 is the 2 to 2 case, channel 3 is the
2 to 3 case, channel 4 is the 2 to 4 case.
"""

import math
from collections import Counter

import numpy as np

DEFAULT_FILE = "events.txt"

# Slot index 0 is "particle 3", because particles 1 and 2 are the incoming beam
# protons and are not written to the file.
FIRST_PARTICLE_NUMBER = 3


class Channel:
    """Every logged event that produced the same number of particles.

    All the arrays have shape (n_events, n_particles). Column j is always the
    j-th particle row as written in the file, so column 0 is particle 3.
    """

    def __init__(self, n_particles, rows, states):
        self.n_particles = n_particles
        self.states = states                      # one name tuple per event

        block = np.asarray(rows, dtype=float).reshape(-1, n_particles, 6)
        self.energy = block[:, :, 0]
        self.theta = block[:, :, 1]               # degrees, 0 to 180
        self.phi = block[:, :, 2]                 # degrees, 0 to 360
        self.px = block[:, :, 3]
        self.py = block[:, :, 4]
        self.pz = block[:, :, 5]

    @property
    def n_events(self):
        return self.energy.shape[0]

    def particle_number(self, slot):
        """File slot index to the particle number used in the report (3, 4, 5, 6)."""
        return slot + FIRST_PARTICLE_NUMBER

    def names_in_slot(self, slot):
        """Every particle name that has appeared in this slot, most common first.

        A slot is not always the same particle. The 2 to 4 channel mixes several
        final states, so slot 0 is sometimes a muon and sometimes a photon. That
        does not change the kinematics, because the generator treats every
        outgoing particle as massless, but it does change what a sensible label
        for the curve is.
        """
        counts = Counter(state[slot] for state in self.states)
        return [name for name, _ in counts.most_common()]

    def label(self, slot):
        """A legend label like 'Particle 3 (photon)' or 'Particle 3 (mixed)'."""
        names = self.names_in_slot(slot)
        tag = names[0] if len(names) == 1 else "mixed"
        return f"Particle {self.particle_number(slot)} ({tag})"

    def four_vector(self, slot):
        """(E, px, py, pz) arrays for one particle slot, one entry per event."""
        return (
            self.energy[:, slot],
            self.px[:, slot],
            self.py[:, slot],
            self.pz[:, slot],
        )

    def invariant_mass(self, *slots):
        """Invariant mass of a group of particles, one value per event.

        Adds the four vectors of the chosen slots together and takes the mass of
        the total. Clamped at zero first, because floating point noise can push
        a genuinely massless combination a hair below zero inside the root. This
        is the 3D form, so p_y is included -- the old 2D scripts left it out.
        """
        E = np.zeros(self.n_events)
        px = np.zeros(self.n_events)
        py = np.zeros(self.n_events)
        pz = np.zeros(self.n_events)
        for slot in slots:
            E += self.energy[:, slot]
            px += self.px[:, slot]
            py += self.py[:, slot]
            pz += self.pz[:, slot]
        m2 = E ** 2 - px ** 2 - py ** 2 - pz ** 2
        return np.sqrt(np.clip(m2, 0.0, None))

    def filter_state(self, state):
        """A new Channel holding only the events with this exact final state.

        Pass a tuple of names, for example ("photon", "photon", "proton",
        "proton"). Useful for pulling the diphoton events out of channel 4.
        """
        state = tuple(state)
        keep = [i for i, s in enumerate(self.states) if s == state]
        if not keep:
            raise ValueError(f"no events in channel {self.n_particles} with state {state}")

        block = np.stack(
            [self.energy, self.theta, self.phi, self.px, self.py, self.pz], axis=2
        )[keep]
        rows = block.reshape(-1, 6).tolist()
        return Channel(self.n_particles, rows, [self.states[i] for i in keep])


def load_events(filename=DEFAULT_FILE):
    """Read the whole events file. Returns {n_particles: Channel}.

    Anything that is not a particle row is skipped by trying to read its six
    numbers and moving on if that fails, so the column header line is filtered
    out without having to match its exact text.
    """
    rows = {}
    states = {}

    n_expected = 0
    event_rows = []
    event_names = []

    def flush():
        # Only keep an event whose row count matches its own N_Particles line,
        # so a half written final event at the end of the file cannot corrupt
        # the reshape further down.
        if n_expected and len(event_rows) == n_expected:
            rows.setdefault(n_expected, []).extend(event_rows)
            states.setdefault(n_expected, []).append(tuple(event_names))

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("___EVENT_ID"):
                flush()
                n_expected = 0
                event_rows = []
                event_names = []
                continue

            if line.startswith("N_Particles"):
                n_expected = int(line.split(":")[1])
                continue

            parts = line.split()
            if len(parts) != 7:
                continue
            try:
                values = [float(v) for v in parts[1:]]
            except ValueError:
                continue  # the "Particle Energy(TeV) ..." column header

            event_names.append(parts[0])
            event_rows.append(values)

    flush()

    return {n: Channel(n, rows[n], states[n]) for n in sorted(rows)}


# The name the report uses for each channel, for titles.
CHANNEL_NAMES = {2: "2 to 2", 3: "2 to 3", 4: "2 to 4"}


def channel_name(n_particles):
    return CHANNEL_NAMES.get(n_particles, f"2 to {n_particles}")


if __name__ == "__main__":
    # Quick look at what is in the file.
    channels = load_events()
    for n, ch in channels.items():
        print(f"{channel_name(n)}: {ch.n_events} events")
        for slot in range(n):
            names = ", ".join(ch.names_in_slot(slot))
            mean_E = ch.energy[:, slot].mean()
            print(f"    slot {slot} (particle {ch.particle_number(slot)}): "
                  f"mean E = {mean_E:7.4f} TeV   names = {names}")
