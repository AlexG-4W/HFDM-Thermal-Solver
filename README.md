# HFDM Thermal Solver v1.1

![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)
![Status](https://img.shields.io/badge/status-Stable-brightgreen.svg)

**HFDM (Heterogeneous Finite Difference Method) Thermal Solver** is a high-performance 2D computational heat transfer simulation tool specifically designed for Printed Circuit Boards (PCBs). 

Unlike standard FDM solvers that assume uniform material properties, HFDM v1.1 parses actual manufacturing Gerber files (RS-274X) to build a heterogeneous thermal conductivity matrix, accurately simulating the thermodynamic behavior at the interfaces between FR-4 substrate and copper traces.

---

## Key Features

* **Direct Gerber Parsing (RS-274X):** Automates the extraction of top-layer copper topology using pcb-tools, mapping it directly into the spatial thermal conductivity matrix K(x,y).
* **Heterogeneous Interface Physics:** Implements Harmonic Mean Averaging on staggered grids. This strictly enforces the conservation of energy and prevents non-physical heat diffusion at the sharp boundaries between high-conductivity copper (385 W/mK) and low-conductivity FR-4 (0.3 W/mK).
* **Strict CFL Stability:** The Transient solver automatically calculates a dynamic time step (dt) based on the Courant-Friedrichs-Lewy condition, utilizing the maximum thermal diffusivity (alpha_max) and convective heat sinks to guarantee absolute mathematical stability (zero numerical explosions).
* **Zero For-Loops (100% Vectorized):** The mathematical core relies entirely on NumPy array slicing. Convective boundary conditions are implemented via the Ghost Node method, preserving O(dx^2) spatial accuracy without iterative bottlenecks.
* **Multithreaded PyQt6 GUI:** A robust desktop interface that isolates the intensive FDM calculations into background QRunnable threads, preventing UI freezing while utilizing pyqtSignal for safe memory management and Matplotlib rendering.

---


<img width="2545" height="1376" alt="scr1" src="https://github.com/user-attachments/assets/7e7fe15d-90be-4f16-97f5-fe19bc6feccc" />

<img width="2559" height="1394" alt="scr2" src="https://github.com/user-attachments/assets/cc3e082d-5ec8-4afd-a6a6-91592efae555" />





## Mathematical Foundation

The core solver is built upon the divergent form of the 2D Heat Equation, accounting for spatially varying thermal conductivity K(x,y), volumetric heat generation Q(x,y), and convective cooling h:

$$ \rho c_p \frac{\partial u}{\partial t} = \nabla \cdot (K(x,y) \nabla u) + Q(x,y) - \frac{2h(u - T_{amb})}{d} $$

To accurately calculate the heat flux between two adjacent cells (e.g., node $i$ and $i+1$) with vastly different materials, the solver computes the effective interface conductivity using the harmonic mean:

$$ K_{interface} = \frac{2 K_i K_{i+1}}{K_i + K_{i+1}} $$

---

## Project Architecture

* `solver.py` — The mathematical engine. Contains the highly vectorized `solve_steady_state` (Jacobi iteration) and `solve_transient` (Explicit Euler) algorithms.
* `data_loader.py` — Handles the parsing of CSV power profiles, Heatsink definitions, and rasterizes Gerber files to build the baseline `k_eff` matrix.
* `visualization.py` — Matplotlib integration for rendering real-time thermal maps and hardware topology.
* `gui.py` — The PyQt6 user interface, signal routing, and application state management.
* `config.py` — Global physical constants, spatial resolution (`dx`), and PCB stackup parameters.

---

## Installation & Usage

1. Clone the repository:
```bash
git clone https://github.com/AlexG-4W/HFDM-Thermal-Solver
cd HFDM-Thermal-Solver
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```
*(Required packages: `numpy`, `matplotlib`, `PyQt6`, `pcb-tools`, `pycairo`)*

3. Run the application:
```bash
python main.py
```

### Running a Test Simulation
Check the `/examples` folder for sample files. 
1. Load `bga_components.csv` to import the heat sources.
2. Load `bga_heatsinks.csv` to apply custom convective zones.
3. Load `bga_top_copper.gbr` to generate the correct thermal pathways.
4. Click **Run Steady-State** to visualize the thermal equilibrium.

---

## License
Distributed under the GNU General Public License v3.0
