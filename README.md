# LHC Subatomic Particle Collision Simulator

A small Monte Carlo simulation of proton–proton collisions at the Large Hadron
Collider (LHC), written in Python. It "collides" two protons over and over,
invents a physically valid set of outgoing particles for each collision, filters
them through a simplified detector, and saves the survivors to a text file that
can then be turned into histograms.

This README is written for someone who has **never seen the code or the physics
before**. It starts from the physics ideas, then shows how each idea turns into
code, and finally how to run everything.

---

## Table of contents

1. [The physics, from scratch](#1-the-physics-from-scratch)
2. [How the code mirrors the physics](#2-how-the-code-mirrors-the-physics)
3. [How to run it](#3-how-to-run-it)
4. [The output file explained](#4-the-output-file-explained)
5. [Making the plots](#5-making-the-plots)
6. [Known limitations / things still to fix](#6-known-limitations--things-still-to-fix)
7. [File map](#7-file-map)

---

## 1. The physics, from scratch

Read this section first. If you understand it, the code in section 2 will feel
like a straight translation of these ideas.

### 1.1 What is actually being simulated?

The LHC is a 27 km ring that accelerates two beams of **protons** in opposite
directions and smashes them head-on. Each beam carries **6.8 TeV** of energy, so
a head-on collision has a total energy of **13.6 TeV** (this is the real LHC
"Run 3" energy).

- **TeV** = tera-electronvolt, a unit of energy. Everything in this project is
  measured in TeV.
- A **proton** is not fundamental — it is a bag of smaller pieces (quarks and
  gluons, collectively "partons"). When two protons collide, it is really one
  piece from each proton that interacts hard, while the leftovers carry on.

When the collision happens, the energy is briefly concentrated into a tiny point
and then re-materialises as **new particles**. Which particles come out is random
and governed by probabilities.

### 1.2 What comes out of a collision? (the "final state")

The list of particles produced by one collision is called the **final state**.
This simulation uses a fixed menu of possible final states, each with a fixed
probability of occurring:

| Probability | Final state (particles produced)              |
|-------------|-----------------------------------------------|
| 25%         | photon + proton                               |
| 18%         | neutron + antineutron + proton + proton       |
| 22%         | photon + proton + proton                      |
| 15%         | positron + electron + proton + proton         |
| 15%         | muon + antimuon + proton + proton             |
| 5%          | photon + photon + proton + proton             |

Most cases also produce "proton remnants" (the leftovers of the original protons
that didn't interact), alongside the interesting particles from the hard
collision (photons, an electron pair, etc.).

The most physically interesting case is the **two-photon** one
(`photon + photon + ...`): two photons are how a **Higgs boson** reveals itself.
It is deliberately rare here (5%), just as a real Higgs signal is rare.

### 1.3 The two rules nature never breaks

Whatever comes out, two quantities must be **conserved** — identical before and
after the collision:

1. **Energy.** The outgoing particles' energies must add up to the total we
   started with, **13.6 TeV**.
2. **Momentum.** Momentum is "quantity of motion" and has a direction. The two
   protons come in exactly head-on with equal and opposite momentum, so the
   total momentum before the collision is **zero**. Therefore the outgoing
   particles' momenta must also add up to **zero** (they fly out in balanced
   directions, like fragments of an explosion).

Every event this simulation produces obeys both rules exactly. That is the
single most important correctness property of the whole program.

### 1.4 Describing one particle: the "four-vector"

To keep the math simple, this simulation lives in a **2D plane** with two
directions: **`z` is the horizontal (left–right) axis** — along the beam — and
**`x` is the vertical (up–down) axis** — sideways to the beam. The emission angle
is measured from the `z`-axis, so 0° points right along the beam and 90° points
straight up.

Each particle is described by three numbers bundled together, called a
**four-vector**:

```
(E, px, pz)   =   (energy, momentum sideways, momentum along the beam)
```

Two derived quantities matter:

- **Momentum magnitude** `|p| = sqrt(px² + pz²)` — how much motion it has,
  ignoring direction.
- **Emission angle** — which direction it flew, computed from `px` and `pz`.

### 1.5 Mass, and the most important formula: invariant mass

Einstein's relation ties a particle's energy, momentum, and mass together:

```
E² = |p|² + m²          (in units where the speed of light = 1)
```

So a particle's **mass** can be recovered from its four-vector:

```
m = sqrt(E² − px² − pz²)
```

This is called the **invariant mass**. "Invariant" means every observer agrees
on it, no matter how fast they are moving — which makes it the perfect tool for
identifying particles.

The magic trick: you can compute the invariant mass of **several particles added
together**, treating them as one object:

```
m(of a group) = sqrt( (ΣE)² − (Σpx)² − (Σpz)² )
```

**Why we care:** if two photons came from a decaying Higgs boson, then the
invariant mass of those two photons will always equal the **Higgs mass, 125 GeV
(0.125 TeV)**. So if you make a histogram of the two-photon invariant mass over
many events, a Higgs shows up as a **bump at 125 GeV** sitting on top of a smooth
background. (This is exactly how the Higgs was discovered in 2012.)
*Note: this version of the simulation does not inject a Higgs bump — see
[section 6](#6-known-limitations--things-still-to-fix).*

### 1.6 How we build one event: split two at a time

Here is the clever idea (from the project report, section 2.5.3) that lets us
create any final state while **guaranteeing** energy and momentum conservation.

Rather than trying to place 3 or 4 particles at once, we only ever split **one
thing into two**, because a two-body split is easy to make conservation-perfect:

> If a parent particle sits still and splits into two, the two children must fly
> off in **exactly opposite directions** with **equal and opposite momentum**,
> and energy conservation fixes exactly how much energy each child gets.

To build a 4-particle final state (say two photons + two protons):

1. Pretend the two photons are secretly **one** made-up particle, call it "(34)",
   and the two protons are another made-up particle, "(56)". Each made-up
   particle has an invariant mass we pick at random.
2. Split the whole collision into `(34)` + `(56)` — one clean two-body split.
3. Now split `(34)` into its two photons, and `(56)` into its two protons — two
   more two-body splits.

There is one complication. When we split `(34)` into two photons, we do it in the
frame where `(34)` is standing still. But `(34)` is actually **flying through the
lab**. So we have to "**boost**" the two photons — a Lorentz transformation that
accounts for that motion, speeding them up in the right direction. This is
standard special relativity (report section 2.5.4).

A 3-particle final state works the same way with one made-up composite; a
2-particle final state is a single split with no boosting needed.

The upshot: every particle comes out **on-shell** (its `E`, `p`, and `m` are
consistent), and the whole event conserves energy and momentum automatically.

### 1.7 The detector: what we can actually see

A real detector cannot see everything. This simulation models two blind spots
with **cuts** (acceptance rules). An event is only kept if **every** particle in
it satisfies both:

1. **Energy cut:** energy ≥ **0.02 TeV** (20 GeV). Too-soft particles are lost in
   the noise.
2. **Angle cut:** the particle's angle must fall in **10°–170°** or **190°–350°**.
   The excluded slivers around 0°/180°/360° are the **beam pipe** direction,
   where real detectors have no coverage.

Events that fail are thrown away; only accepted events are saved.

---

## 2. How the code mirrors the physics

The project is split into **two Python files** that map directly onto the two
halves of section 1:

| File             | Plays the role of… | Responsible for |
|------------------|--------------------|-----------------|
| `kinematics.py`  | the **physicist**  | the physics: four-vectors, two-body splits, boosts, invariant mass (section 1.4–1.6). Deals only in numbers, never particle names. |
| `simulation.py`  | the **director**   | running the experiment: choosing final states, applying detector cuts, and writing the output (section 1.2, 1.3, 1.7). |

### 2.1 `kinematics.py` — the physics engine

It is organised in five numbered sections, in reading order:

1. **The four-vector.** `P4 = (E, px, pz)` and the functions `mass(p4)` and
   `invariant_mass(*p4s)` — a direct translation of section 1.4 and 1.5.
2. **The two building blocks.**
   - `two_body_decay(parent_mass, m1, m2)` — splits a parent-at-rest into two
     children (section 1.6, the "one thing into two" rule).
   - `boost(p4, parent_p4)` — the Lorentz boost that carries a child from its
     parent's rest frame into the lab.
   - `decay_in_lab(...)` — a tiny helper that does "decay, then boost both
     children," so the generators below don't repeat themselves.
3. **`draw_composite_mass(low, high)`** — randomly picks the invariant mass of a
   made-up composite like `(34)`.
4. **The generators** — one per final-state size, named after the report's
   cases: `two_to_two`, `two_to_three`, `two_to_four`. Each is just the building
   blocks arranged as described in section 1.6. All the *final* particles are
   treated as **massless**; only the intermediate composites carry mass.
5. **`generate_momenta(n_particles, total_energy)`** — the single entry point.
   Since every particle is massless, it only needs to know *how many* particles
   there are; it calls the right generator and returns one `P4` per particle.

`kinematics.py` never mentions "photon" or "proton" — it only sees masses and
numbers, which keeps the physics reusable and testable on its own.

### 2.2 `simulation.py` — the experiment

It reads top-to-bottom as five steps:

1. **`choose_final_state()`** — rolls a random number and returns the list of
   particle **names** for this collision (the menu from section 1.2).
2. **`make_event(names)`** — the bridge to the physics. It calls
   `kinematics.generate_momenta(...)` (passing just how many particles there are,
   since they're all massless) and packages the results into a list of
   **`Particle`** objects. A `Particle` bundles everything about one outgoing
   particle: `(name, energy, px, pz, angle)`.
3. **The detector** — three small functions:
   - `is_seen(particle)` — the energy + angle cuts for one particle (section 1.7).
   - `passes_cuts(event)` — keeps the event only if **every** particle is seen.
   - `is_conserved(event)` — a safety re-check that energy sums to 13.6 and
     momentum sums to zero (section 1.3).
4. **`write_event(...)`** — writes one accepted event to the output file.
5. **`run_simulation()` / `main()`** — the loop: keep colliding until 100 good
   events have been logged, then print a short summary.

### 2.3 The data flow, end to end

```
 simulation.py                                   kinematics.py
 ─────────────                                   ─────────────
 choose_final_state()
   → ["photon","photon","proton","proton"]   (just names, no physics yet)
        │
   count them  → 4 particles
        │
   make_event()  ── generate_momenta(4, 13.6) ──►  pick generator (two_to_four)
                                                     split into composites,
                                                     decay each, boost
        ◄──────────────  [P4, P4, P4, P4]  ──────────────
        │
   wrap each P4 + name + angle  →  [Particle, Particle, Particle, Particle]
        │
   is_conserved? ✓   passes_cuts? ✓
        │
   write_event()  →  append to events.txt
        │
   repeat until 100 events are logged
```

**One sentence:** `simulation.py` decides *which* particles appear and judges the
result; `kinematics.py` produces the *actual energies and momenta*, obeying the
conservation laws by construction.

---

## 3. How to run it

### 3.1 What you need

- **Python 3** (this project was developed with the MSYS2 / UCRT64 Python at
  `E:\msys64\ucrt64\bin\python.exe`).
- **matplotlib** — only needed for the plotting scripts, not for the core
  simulation.

> **Installing matplotlib on this setup:** the MSYS2 Python does **not** use
> `pip`. Install packages with its own package manager instead:
> ```
> pacman -S mingw-w64-ucrt-x86_64-python-matplotlib
> ```
> If VS Code cannot find matplotlib even after installing, make sure it is using
> the same interpreter: **Ctrl+Shift+P → "Python: Select Interpreter" →
> `E:\msys64\ucrt64\bin\python.exe`**.

### 3.2 Run the simulation

From the project folder:

```
python simulation.py
```

This runs collisions until **100** events pass the detector, writes them to
`events.txt`, and prints something like:

```
Simulation finished.
Total collisions: 153
Accepted events: 100
Saved to: events.txt
```

("Total collisions" is larger than "Accepted events" because many collisions are
thrown out by the detector cuts.)

You can change the behaviour by editing the settings at the top of
`simulation.py`:

- `TARGET_LOGGED_EVENTS` — how many accepted events to collect (e.g. raise to
  `50000` for smooth histograms).
- `MIN_ENERGY`, `VISIBLE_ANGLES_DEG` — the detector cuts.
- `OUTPUT_FILE` — the output filename.

---

## 4. The output file explained

`events.txt` contains one block per accepted event:

```
___EVENT_ID: 1___
N_Particles: 4
Particle       Energy(TeV)    Angle(deg)     p_x(TeV)       p_z(TeV)
photon         2.537252       39.647730      1.618933       1.953638
photon         5.667125       218.592597     -3.535031      -4.429431
proton         4.019932       48.199762      2.914029       2.605465
proton         1.375691       262.596430     -0.997931      -0.129672

___EVENT_ID: 2___
...
```

- `___EVENT_ID: n___` — event number.
- `N_Particles` — how many particles are in this event (2, 3, or 4).
- Then one row per particle: its **name**, **energy** (TeV), **angle**
  (degrees), and the two momentum components **p_x**, **p_z** (TeV).

You can sanity-check any event by hand: the four energies add up to 13.6, and the
four `p_x` values (and the four `p_z` values) each add up to ~0.

---

## 5. Making the plots

Three analysis scripts read the saved events and draw histograms with matplotlib:

- `graphing.py` — an energy histogram for each particle position.
- `angle_graphing.py` — an angle (cos θ) histogram for each particle position.
- `invariant_mass_graph_1&2.py` — the two-photon **invariant mass** histogram
  (the plot where a Higgs would appear as a bump).

Run one like this (the `&` in the filename must be quoted):

```
python "invariant_mass_graph_1&2.py"
```

> **Important — these scripts are currently out of sync with the simulation.**
> They were written for an older output file named `no_higgs_events.txt` with
> angles in **radians**, but `simulation.py` now writes `events.txt` with angles
> in **degrees**, and `events.txt` mixes 3- and 4-particle events. Before the
> plots will work you must update each script's `DATA_FILE` to `"events.txt"`,
> and make sure it only reads the 4-particle (two-photon) events. See section 6.

---

## 6. Known limitations / things still to fix

These are honest caveats, not bugs that break the run:

1. **Every particle is treated as massless.** This is a deliberate
   simplification: the final particles get no rest mass, so the `MASS` table in
   `simulation.py` is kept only for reference and isn't used when generating
   events. It's a good approximation because all these masses are tiny next to
   13.6 TeV, but it does mean a reconstructed proton comes out with ~0 mass, not
   its real 0.938 GeV. (The intermediate composites still carry mass — that part
   is essential and unaffected.)
2. **The plotting scripts need updating** to read `events.txt` (degrees) instead
   of `no_higgs_events.txt` (radians) — see section 5.
3. **No Higgs bump.** The two-photon invariant mass is drawn as a plain random
   spread, so there is no 125 GeV peak. Adding one back means drawing the `(34)`
   composite mass from a narrow Gaussian at 0.125 TeV for a fraction of events.
4. **Energy sharing is random, not from real physics.** The composite masses are
   drawn uniformly, not from real particle-physics probabilities ("matrix
   elements"), so the distributions are illustrative, not predictive.
5. **2D only.** The simulation works in an (x, z) plane, not full 3D, to keep the
   math approachable.

---

## 7. File map

| File                         | What it is |
|------------------------------|------------|
| `simulation.py`              | Runs the experiment: choose final states, apply detector cuts, write `events.txt`. **Start here / run this.** |
| `kinematics.py`              | The physics engine: four-vectors, two-body decays, Lorentz boosts, invariant mass. |
| `events.txt`                 | The output — one block per accepted event (created when you run the simulation). |
| `graphing.py`                | Energy histograms per particle (needs updating, section 5). |
| `angle_graphing.py`          | Angle histograms per particle (needs updating, section 5). |
| `invariant_mass_graph_1&2.py`| Two-photon invariant-mass histogram (needs updating, section 5). |
| `LHC_Simulation_Report.pdf`  | The written report with the full physics derivations (section 2.5 = the engine). |

---

*Phenikaa University Research Lab — computational physics simulation of
high-energy proton collisions.*
