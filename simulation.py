"""LHC collision simulation.

This file runs the experiment. It decides which particles come out of each
collision, then judges the result with a simple detector and logs the events
that pass. All the actual physics, meaning the energies, momenta, conservation,
and boosts, lives in kinematics.py.

The event flow from top to bottom is choose_final_state, then make_event, then
the detector cuts, then write to file.
"""

import math
import random
from collections import namedtuple

import kinematics

TOTAL_ENERGY = 13.6           # TeV, the LHC collision energy
TARGET_LOGGED_EVENTS = 100    # stop once this many events pass the detector
OUTPUT_FILE = "events.txt"

# Detector cuts. A particle is only seen if it is energetic enough and does not
# disappear down the beam pipe.
MIN_ENERGY = 0.02                             # TeV
VISIBLE_ANGLES_DEG = [(10, 170), (190, 350)]  # angle windows the detector covers

# Real particle rest masses in TeV. Kept for reference only. The simulation now
# treats every outgoing particle as massless, so these are not used when
# generating events. They would matter again if we ever switch back to real masses.
MASS = {
    "photon": 0.0,
    "proton": 0.000938,
    "positron": 0.000000511,
    "electron": 0.000000511,
    "muon": 0.0001057,
    "antimuon": 0.0001057,
    "neutron": 0.000939,
    "antineutron": 0.000939,
}

# One outgoing particle. It carries its name plus everything the detector
# measures about it.
Particle = namedtuple("Particle", ["name", "energy", "px", "pz", "angle"])


# Step 1. Choose what comes out of the collision.
def choose_final_state():
    """Randomly pick the list of particles produced by one collision."""
    r = random.random()
    if r < 0.25:
        return ["photon", "proton"]
    elif r < 0.43:
        return ["neutron", "antineutron", "proton", "proton"]
    elif r < 0.65:
        return ["photon", "proton", "proton"]
    elif r < 0.80:
        return ["positron", "electron", "proton", "proton"]
    elif r < 0.95:
        return ["muon", "antimuon", "proton", "proton"]
    else:
        return ["photon", "photon", "proton", "proton"]


# Step 2. Turn that list of names into a real event.
def make_event(names):
    """Give each named particle its energy, momentum, and angle.

    We ask the physics engine for the four momenta. It only needs to know how
    many particles there are, since every particle is treated as massless. Then
    we pair each result back with its name and work out its emission angle. We
    return a list of Particle objects.
    """
    momenta = kinematics.generate_momenta(len(names), TOTAL_ENERGY)
    if momenta is None:
        return None

    event = []
    for name, p in zip(names, momenta):
        # z is the horizontal left right axis and x is the vertical up down axis.
        # The angle is measured from the z axis, so 0 degrees points right along z
        # and 90 degrees points straight up along x.
        angle = math.atan2(p.px, p.pz) % (2 * math.pi)
        event.append(Particle(name, p.E, p.px, p.pz, angle))
    return event


def collision():
    """Simulate one whole collision. Choose a final state, then generate it."""
    return make_event(choose_final_state())


# Step 3. The detector.
def is_seen(particle):
    """True if the detector can measure this particle, meaning it passes the cuts."""
    if particle.energy < MIN_ENERGY:
        return False
    angle_deg = math.degrees(particle.angle)
    return any(low <= angle_deg <= high for low, high in VISIBLE_ANGLES_DEG)


def passes_cuts(event):
    """The event is kept only if every particle in it is seen."""
    return all(is_seen(p) for p in event)


def is_conserved(event):
    """Sanity check that energy sums to 13.6 TeV and momentum sums to zero.

    The physics engine guarantees this by construction. We re-check to catch any
    mistake. The small tolerance absorbs the floating point noise that Lorentz
    boosts leave behind.
    """
    total_E = sum(p.energy for p in event)
    total_px = sum(p.px for p in event)
    total_pz = sum(p.pz for p in event)
    return (
        math.isclose(total_E, TOTAL_ENERGY, abs_tol=1e-6)
        and math.isclose(total_px, 0.0, abs_tol=1e-6)
        and math.isclose(total_pz, 0.0, abs_tol=1e-6)
    )


# Step 4. Write an accepted event to the file.
def write_event(f, event_id, event):
    f.write(f"___EVENT_ID: {event_id}___\n")
    f.write(f"N_Particles: {len(event)}\n")
    f.write(
        f"{'Particle':<15}{'Energy(TeV)':<15}{'Angle(deg)':<15}"
        f"{'p_x(TeV)':<15}{'p_z(TeV)':<15}\n"
    )
    for p in event:
        f.write(
            f"{p.name:<15}{p.energy:<15.6f}{math.degrees(p.angle):<15.6f}"
            f"{p.px:<15.6f}{p.pz:<15.6f}\n"
        )
    f.write("\n")


# Step 5. Run collisions until enough events pass.
def run_simulation():
    """Collide until TARGET_LOGGED_EVENTS good events are logged. Returns counts."""
    logged = 0
    total = 0
    with open(OUTPUT_FILE, "w") as f:
        while logged < TARGET_LOGGED_EVENTS:
            total += 1
            event = collision()

            # Skip anything impossible, not conserving, or unseen by the detector.
            if event is None or not is_conserved(event) or not passes_cuts(event):
                continue

            logged += 1
            write_event(f, logged, event)
    return total, logged


def main():
    total, logged = run_simulation()
    print("Simulation finished.")
    print(f"Total collisions: {total}")
    print(f"Accepted events: {logged}")
    print(f"Saved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
