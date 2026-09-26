"""Energy of every outgoing particle, drawn as outlines on shared axes.

One figure per channel. Every particle in that channel gets a curve, so you can
see straight away which particle carries more of the 13.6 TeV.

Replaces energy_graphing.py and energy_graphing_2to3.py, which drew one filled
histogram per particle in separate files.

Run:  python graph_energy_steps.py
"""

import os
from pathlib import Path

import event_data
import step_plot

# Resolved from this file's location, not the working directory, so the
# script runs the same from anywhere.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_FILE = PROJECT_ROOT / "collision_data" / "events.txt"
OUTPUT_DIR = PROJECT_ROOT / "plots"

TOTAL_ENERGY = 13.6   # TeV, matches simulation.TOTAL_ENERGY
BINS = 50


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    channels = event_data.load_events(DATA_FILE)

    for n, ch in channels.items():
        name = event_data.channel_name(n)
        print(f"{name} case: {ch.n_events} events")

        series = [(ch.label(slot), ch.energy[:, slot]) for slot in range(n)]

        step_plot.outline_figure(
            series,
            bins=BINS,
            value_range=(0.0, TOTAL_ENERGY),
            xlabel="Energy (TeV)",
            title=f"Energy of each particle, {name} case",
            subtitle=f"{ch.n_events} events, {TOTAL_ENERGY} TeV total",
            filename=os.path.join(OUTPUT_DIR, f"energy_steps_{n}particles.png"),
        )

        # The 2 to 2 case has both particles pinned to exactly half the collision
        # energy, so its two curves are a single spike sitting on top of each
        # other. That is correct, not a drawing bug -- with two massless
        # particles and nothing else, energy conservation leaves no freedom.
        if n == 2:
            print("    note: both curves are one spike at 6.8 TeV, they overlap exactly")


if __name__ == "__main__":
    main()
