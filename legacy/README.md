# legacy/

The graphing scripts from before the simulation moved to 3D, kept for reference.
Nothing in the current pipeline imports any of them.

They are archived rather than deleted because they are the scripts the earlier
results in the report were made with, and because `energy_graphing_2to3_E45.py`
is the one that settled the "why is the E3 histogram curving" question.

## Why they were replaced

Every one of these scripts contains its own copy of the same file parser, and
each copy hard codes a **column number**. When the simulation went to 3D the row
layout gained two columns:

    before   name  Energy  Angle(rad)  p_x  p_z
    after    name  Energy  Theta(deg)  Phi(deg)  p_x  p_y  p_z

Nothing errors when a column number goes stale. The scripts keep running and
quietly plot the wrong quantity. That is the whole reason the current pipeline
reads the file in one place, `event_data.py`, and refers to columns by name.

## Status of each file

| File | Still correct? | Notes |
|---|---|---|
| `energy_graphing.py` | yes | Reads column 1, which is Energy in both layouts |
| `energy_graphing_2to3.py` | yes | Same |
| `energy_graphing_2to3_E45.py` | yes | Same. Plots E45 = 13.6 - E3, the composite energy. Was saved without a `.py` extension as `2to3_energy_particle(45)` |
| `graphing.py` | **no** | Column 2 used to be the planar angle over 0 to 360. It is now Theta, which stops at 180, still plotted over `range=(0, 360)`. Half the axis is empty and phi is never plotted |
| `graphing_2to3.py` | **no** | Same fault |
| `invariant_mass_graph_1&2.py` | only on its own data | Reads `parts[3]` as p_x and `parts[4]` as p_z. In the 3D layout those are **Phi** and **p_x**. It still works because it points at `no_higgs_events.txt`, which is old format. Pointed at `events.txt` it gets a negative m^2 on 81959 of 81962 groups and draws an empty histogram |

## Running them

They resolve `events.txt` relative to the working directory, so run them from the
project root, not from inside this folder:

    python legacy/energy_graphing.py

They also write their PNGs to the working directory using bare filenames, so
running one from the root drops loose images back into the root. The old outputs
are already saved in `legacy/plots/`.

## Current equivalents

| Old | New |
|---|---|
| `energy_graphing.py`, `energy_graphing_2to3.py` | `graph_energy_steps.py` |
| `graphing.py`, `graphing_2to3.py` | `graph_angle_steps.py` |
| `invariant_mass_graph_1&2.py` | `graph_mass_steps.py` |

The new scripts cover all three channels, including 2 to 4, which is the
majority of the events and which none of the old scripts handled.
