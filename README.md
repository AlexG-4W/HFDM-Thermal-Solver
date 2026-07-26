# HFDM Thermal Solver 


A 2D finite-difference thermal simulator for printed circuit boards. It solves
the steady-state and transient heat equation over a board whose thermal
conductivity is derived from the real copper topology, so component placement,
copper pours and heatsink zones can be evaluated before hardware exists.

The board is modelled as a thin plate: heat spreads laterally through the
copper/FR-4 stackup and leaves through both faces by convection and radiation.

---

## Key features

### Physics

- **Convection and radiation.** Radiative loss uses the exact factorisation of
  the Stefan-Boltzmann law, `εσ(T⁴ − Ta⁴) = h_rad·(T − Ta)` with
  `h_rad = εσ(T² + Ta²)(T + Ta)`. Because `h_rad` depends on the field being
  solved for, the steady-state solve iterates until it stops moving. This is
  not a refinement: on a small board at 110 °C radiation carries **56 % of the
  heat load**, and omitting it makes a board look roughly twice as hot as it is.
- **Stackup-aware conductivity.** Every cell is a through-thickness composite,
  never bulk copper: `k_cell = k_eff(layers, oz) + coverage · k_delta`. A
  4-layer 1 oz board spans 33.96–42.38 W/mK.
- **Energy conservation.** Measured closure between generated power and the sum
  of face, edge and radiative losses is **0.064 %** on a real 66 × 35 mm board.
- **Two solvers.** Steady state uses red-black SOR with the boundary strips
  included in the colouring, which allows ω ≈ 1.96 and converges ~12.7× faster
  than the naive scheme. Transient uses explicit Euler on a CFL-limited step,
  with a divergence guard.
- **Honest convergence reporting.** The solver returns a convergence flag; a run
  that hits the iteration cap says so instead of presenting a partial answer as
  a result.


<img width="1600" height="950" alt="1" src="https://github.com/user-attachments/assets/4f96c61b-3718-462d-820c-c025c651f08a" />


<img width="1600" height="950" alt="2" src="https://github.com/user-attachments/assets/330da65d-7ae6-488c-b0b8-ee5423ed2b9d" />


<img width="1600" height="950" alt="3" src="https://github.com/user-attachments/assets/73f8541d-6e7e-4bc6-a61e-c90e9765b6ef" />





  

### Gerber import

- **X2/X3 via PyGerber**, millimetre and inch sources alike.
- **Geometric registration.** The raster is placed at its true millimetre
  coordinates rather than stretched to fill the board, so copper lines up with
  the component coordinates from the CSV.
- **BOX (area-average) resampling.** Downsampling ~40 px/mm to a 0.5 mm grid
  with nearest-neighbour deletes whole traces; area-averaging preserves them and
  yields a *fractional* copper coverage per cell, so a pour reads denser than a
  thin track.
- **Board frame with an origin.** CAD exports do not place the board at the
  sheet origin — an EasyEDA or Altium board commonly lands at negative Y. The
  frame has an origin as well as a size, and is auto-fitted to the loaded data.

### Visualisation

The canvas is a layer system: what is drawn is a function of state, not of the
last button pressed.

| Base field (exclusive) | Overlays (independent) |
|---|---|
| Temperature | Component outlines |
| Copper conductivity | Copper mask (translucent, coverage-weighted) |
| Power density | Heatsink zones (contours) |
| None | Probes |

- **Interactive virtual probes.** Left-click to place (up to 10), right-click to
  remove the nearest. Probes read the field currently displayed and track a
  transient live.
- **Colour-scale lock** so the same colour means the same temperature on every
  frame of an animation.
- **PNG/JPEG export** of exactly what is on screen, layers included.

### Input robustness

- CSV files from Excel parse regardless of BOM, delimiter (`,` `;` tab), or
  cp1251/UTF-8 encoding. Errors name the file, the line and the column.
- Components that fall outside the board are reported by name with the power
  they take with them, instead of vanishing silently.
- Component power is normalised over the discretised footprint, so the total
  power on the grid equals the nameplate total exactly.

---

## Installation

```bash
git clone https://github.com/AlexG-4W/HFDM-Thermal-Solver.git
cd HFDM-Thermal-Solver
pip install -r requirements.txt
```

Requires **Python 3.9+** (developed and tested on 3.14).

`requirements.txt` pins two dependencies deliberately:

| Pin | Why |
|---|---|
| `pygerber>=2.4,<3` | The importer uses the `pygerber.gerberx3.api` surface (`Rasterized2DLayer`, `LayerProperties.gerber_bounding_box`), which is 2.x-only. |
| `Pillow>=9.1.0` | `Image.Resampling` — the BOX resampling the copper mask depends on — was introduced in 9.1. |
| `matplotlib>=3.8.0` | `ContourSet` became an Artist in 3.8; the heatsink layer calls `remove()` and `set_visible()` on it. |

Run the application:

```bash
python main.py
```

Run the test suite:

```bash
python -m pytest tests -q
```

---

## Usage

1. **Load Components CSV** — heat sources. The board frame is fitted around them
   automatically; the console reports the frame it chose.
2. **Load Heatsinks CSV** *(optional)* — local zones of enhanced convection.
   Independent of the component layer; loading one never discards the other.
3. **Load Top Copper (Gerber)** *(optional)* — parsed on a background thread.
   The frame is re-fitted if the copper extends past it.
4. **Set the environment** in *Parameters*: ambient temperature, convection `h`,
   emissivity `ε`, and — if the automatic frame is not what you want — the board
   origin and size. `Fit to data` recomputes the frame at any time.
5. **Run Steady-State** or **Run Transient**. Both run off the GUI thread and can
   be stopped; closing the window shuts the worker down cleanly.
6. **Inspect.** Toggle layers under the canvas, click the board to place probes,
   and use *Save Result Image* to export.

### Choosing `h` and `ε`

The single biggest lever on the answer is the surface coefficient, not the mesh.

| Condition | `h` [W/m²K] |
|---|---|
| Still air | 5–12 |
| Gentle airflow (~1 m/s) | 20–30 |
| Forced air (2–3 m/s) | 40–80 |

Set `ε = 0.90` for solder mask or bare FR-4; `ε = 0` disables radiation and
reproduces a convection-only model. Heatsink zones from the CSV override the
global `h` locally.

---

## CSV formats

Coordinates are in millimetres and may be negative — they are interpreted in the
CAD sheet frame, the same one the Gerber uses. `Length_mm` is accepted under its
legacy name `Height_mm`.

### Components (`components_power.csv`)

| Designator | Center_X_mm | Center_Y_mm | Width_mm | Length_mm | Power_Watts |
|---|---|---|---|---|---|
| U1_BGA | 50.0 | 50.0 | 35.0 | 35.0 | 8.0 |
| U2_DDR | 75.0 | 70.0 | 12.0 | 8.0 | 1.5 |

### Heatsinks (`heatsinks.csv`)

| Designator | Center_X_mm | Center_Y_mm | Width_mm | Length_mm | Convection_H |
|---|---|---|---|---|---|
| HS1_BGA | 50.0 | 50.0 | 40.0 | 40.0 | 120.0 |

---

## Project layout

| File | Role |
|---|---|
| `main.py` | Entry point. |
| `gui.py` | Window, workers, data ownership. Does not draw. |
| `board_view.py` | Owns the Matplotlib figure and every artist; layer registry and redraw lifecycle. |
| `solver.py` | `PCBSolver`: steady-state SOR, transient Euler, radiation. |
| `data_loader.py` | CSV parsing, power matrix, Gerber rasterisation and projection. |
| `config.py` | Materials, stackup, board frame, environment defaults. |
| `visualization.py` | Standalone plotting helpers (Qt-free). |
| `tests/` | Test suite plus the BGA and SH1 board fixtures. |

---

## Model and its limits

Worth knowing before trusting a number:

- **2D shell.** Through-thickness gradients are not resolved; the board is a
  plate with an effective in-plane conductivity. Reasonable for 1–2 mm boards,
  not for thick metal-core substrates.
- **Component power is a surface source**, spread uniformly over the footprint.
  Package thermal resistance, die size and junction temperature are not modelled
  — the result is board temperature under the part, not `T_j`.
- **Uniform ambient.** No air flow field, no board-to-board or enclosure
  coupling, no view factors — radiation is to an ambient at `T_amb`.
- **Conductivity comes from one copper layer.** Inner planes are folded into
  `k_eff` uniformly; a plane split under a hot part is not seen.

The dominant uncertainty in practice is the input, not the solver: `h` and the
per-component power budget move the answer far more than the mesh does.

---

