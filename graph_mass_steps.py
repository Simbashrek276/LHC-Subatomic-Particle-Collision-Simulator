"""Invariant masses, drawn as outlines on shared axes.

A word on what "the mass of each particle" means here, because it is not what
you might expect. kinematics.py treats every outgoing particle as massless. So
the mass of a single particle is zero by construction, in all three channels,
always. Plotting it tells you nothing about the physics -- but it is still worth
one figure, because it tells you something about the data file: whatever width
that spike has is pure numerical noise, and it sets the floor on how sharp any
other mass peak in this file can possibly be.

The masses that actually carry information are the ones belonging to groups of
particles. Two massless photons flying apart have a real, heavy combined mass.
That is what the composites in kinematics.py are, and it is what a real detector
reconstructs when it looks for a Higgs. So the rest of this script plots the
invariant mass of every pair.

Which pairs matter depends on the channel.

    2 to 2   only one pair exists and it is the whole collision, so its mass is
             pinned at 13.6 TeV. Included for completeness.
    2 to 3   particles 4 and 5 came out of a composite, so m(45) is the mass the
             generator drew. m(34) and m(35) are combinations that were never a
             real object, and they are the "background" shape for comparison.
    2 to 4   particles 3 and 4 came out of one composite and 5 and 6 out of the
             other, so m(34) and m(56) are the drawn masses. The four crossed
             pairs get their own figure.

This also replaces "invariant_mass_graph_1&2.py", which is written against the
old 2D file layout: it reads column 4 as p_z when in the 3D file column 4 is
p_y, and it never picks up p_z at all. Pointed at the current events.txt it
returns wrong masses rather than failing.

Run:  python graph_mass_steps.py
"""

import os
from itertools import combinations

import event_data
import step_plot

DATA_FILE = "events.txt"
OUTPUT_DIR = "plots"

TOTAL_ENERGY = 13.6   # TeV
BINS = 50

# Range for the single particle mass figure. Every value in here should be
# numerical noise, so the scale is small: 0.005 TeV is 5 GeV. That is where the
# spread actually lands, because events.txt stores six decimal places and a
# rounding of 1e-6 TeV on E and p works out to a few GeV once it goes through
# m = sqrt(E^2 - p^2) at 13.6 TeV.
NOISE_RANGE = (0.0, 0.005)


def pair_label(ch, a, b):
    return f"Pair {ch.particle_number(a)}{ch.particle_number(b)}"


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    channels = event_data.load_events(DATA_FILE)

    for n, ch in channels.items():
        name = event_data.channel_name(n)
        print(f"{name} case: {ch.n_events} events")

        # 1. Single particle masses. Expected to be zero everywhere.
        print("  single particle mass (expected zero, so this measures file precision):")
        step_plot.outline_figure(
            [(ch.label(slot), ch.invariant_mass(slot)) for slot in range(n)],
            bins=BINS,
            value_range=NOISE_RANGE,
            xlabel="Reconstructed single particle mass (TeV)",
            title=f"Single particle mass, {name} case",
            subtitle="every outgoing particle is massless by construction, "
                     "so this width is rounding in events.txt",
            filename=os.path.join(OUTPUT_DIR, f"mass_single_steps_{n}particles.png"),
            ratio_panel=False,
        )

        # 2. Pair masses. Split into two figures for the 2 to 4 case, because
        #    six pairs is more curves than can be told apart on one axes.
        pairs = list(combinations(range(n), 2))

        if n == 4:
            # The two pairs that were really composites, then the four that
            # were not. Keeping them apart is the whole point of the figure.
            groups = [
                ("composite pairs", [(0, 1), (2, 3)], "mass_pairs_composites",
                 "pairs 34 and 56 each came out of one composite, so these are "
                 "the masses the generator drew"),
                ("crossed pairs", [(0, 2), (0, 3), (1, 2), (1, 3)], "mass_pairs_crossed",
                 "these pairs were never a single object, so this is the "
                 "combinatorial background shape"),
            ]
        else:
            groups = [("every pair", pairs, "mass_pairs", f"{ch.n_events} events")]

        for group_name, group_pairs, stem, note in groups:
            print(f"  {group_name}:")
            step_plot.outline_figure(
                [(pair_label(ch, a, b), ch.invariant_mass(a, b)) for a, b in group_pairs],
                bins=BINS,
                value_range=(0.0, TOTAL_ENERGY),
                xlabel="Invariant mass of the pair (TeV)",
                title=f"Pair invariant mass, {group_name}, {name} case",
                subtitle=note,
                filename=os.path.join(OUTPUT_DIR, f"{stem}_steps_{n}particles.png"),
            )

    # 3. The diphoton mass, on its own, from the events that really have two
    #    photons. This is the search channel the old script was aiming at.
    if 4 in channels:
        diphoton_state = ("photon", "photon", "proton", "proton")
        try:
            sub = channels[4].filter_state(diphoton_state)
        except ValueError as exc:
            print(f"diphoton: skipped, {exc}")
            return

        print(f"diphoton events: {sub.n_events}")
        step_plot.outline_figure(
            [("Pair 34 (the two photons)", sub.invariant_mass(0, 1)),
             ("Pair 56 (the two protons)", sub.invariant_mass(2, 3))],
            bins=BINS,
            value_range=(0.0, TOTAL_ENERGY),
            xlabel="Invariant mass (TeV)",
            title="Diphoton invariant mass",
            subtitle=f"{sub.n_events} photon photon proton proton events",
            filename=os.path.join(OUTPUT_DIR, "mass_diphoton_steps.png"),
        )


if __name__ == "__main__":
    main()
