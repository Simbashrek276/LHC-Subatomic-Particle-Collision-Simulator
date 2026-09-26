"""Emission angles of every outgoing particle, drawn as outlines on shared axes.

In 3D a particle has two angles, so this writes two figures per channel.

    theta  the polar angle away from the beam line (z axis), 0 to 180 degrees
    phi    the azimuthal angle around the beam, in the x-y plane, 0 to 360

The two do not look the same and are not supposed to.

phi should come out flat. There is no preferred direction around the beam, so
every value is equally likely.

theta should NOT come out flat, and this is the part that changed when the
simulation went from 2D to 3D. Directions are drawn evenly over the surface of a
sphere, and a sphere has more surface area around its equator than near its
poles. A band of angles near theta = 90 is simply bigger than the same band near
theta = 0, so more particles land in it. The shape you should see is a sin(theta)
arch peaking at 90 degrees. A flat theta histogram would mean the sampling is
wrong.

Because of that, this also writes a cos(theta) figure. cos(theta) undoes the
geometry exactly, so that one IS expected to be flat, which makes it the easier
plot to check the sampling against by eye.

Replaces graphing.py and graphing_2to3.py, which read the angle column and drew
it over a 0 to 360 range. In the 3D file that column is theta, which only ever
reaches 180, so those plots have an empty right half and never show phi at all.

Run:  python graph_angle_steps.py
"""

import os
from pathlib import Path

import numpy as np

import event_data
import step_plot

# Resolved from this file's location, not the working directory, so the
# script runs the same from anywhere.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_FILE = PROJECT_ROOT / "collision_data" / "events.txt"
OUTPUT_DIR = PROJECT_ROOT / "plots"

THETA_BINS = 45    # 4 degrees per bin over 0 to 180
PHI_BINS = 45      # 8 degrees per bin over 0 to 360
COS_BINS = 40


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    channels = event_data.load_events(DATA_FILE)

    for n, ch in channels.items():
        name = event_data.channel_name(n)
        print(f"{name} case: {ch.n_events} events")

        # Polar angle from the beam axis.
        print("  theta (polar angle from the beam):")
        step_plot.outline_figure(
            [(ch.label(slot), ch.theta[:, slot]) for slot in range(n)],
            bins=THETA_BINS,
            value_range=(0.0, 180.0),
            xlabel="Polar angle theta from the beam axis (degrees)",
            title=f"Polar angle of each particle, {name} case",
            subtitle=f"{ch.n_events} events, expect a sin(theta) arch peaking at 90 degrees",
            filename=os.path.join(OUTPUT_DIR, f"theta_steps_{n}particles.png"),
        )

        # Azimuthal angle around the beam axis.
        print("  phi (azimuthal angle around the beam):")
        step_plot.outline_figure(
            [(ch.label(slot), ch.phi[:, slot]) for slot in range(n)],
            bins=PHI_BINS,
            value_range=(0.0, 360.0),
            xlabel="Azimuthal angle phi around the beam axis (degrees)",
            title=f"Azimuthal angle of each particle, {name} case",
            subtitle=f"{ch.n_events} events, expect a flat distribution",
            filename=os.path.join(OUTPUT_DIR, f"phi_steps_{n}particles.png"),
        )

        # cos(theta), the flat-if-correct version of the theta plot.
        print("  cos(theta) (flat if the sphere sampling is right):")
        step_plot.outline_figure(
            [(ch.label(slot), np.cos(np.radians(ch.theta[:, slot]))) for slot in range(n)],
            bins=COS_BINS,
            value_range=(-1.0, 1.0),
            xlabel="cos(theta)",
            title=f"cos of the polar angle, {name} case",
            subtitle=f"{ch.n_events} events, expect a flat distribution if sampling is isotropic",
            filename=os.path.join(OUTPUT_DIR, f"costheta_steps_{n}particles.png"),
        )


if __name__ == "__main__":
    main()
