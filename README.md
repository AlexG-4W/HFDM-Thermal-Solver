# HFDM Thermal Solver (Heterogeneous Finite Difference Method)

<img width="1781" height="1233" alt="image" src="https://github.com/user-attachments/assets/49589ced-e1c6-43e8-af23-ae64d79ea5f4" />
<img width="1784" height="1251" alt="image" src="https://github.com/user-attachments/assets/f14aca52-9951-498f-9cfb-5f6b4bbf69a7" />



A high-performance, 2D Finite Difference Method (FDM) thermal solver designed specifically for Printed Circuit Boards (PCBs). This tool simulates both **transient** (dynamic heating/cooling) and **steady-state** thermal distribution across heterogeneous materials, allowing engineers to identify thermal bottlenecks and optimize component placement and copper pours.


## 🚀 Key Features

* **Heterogeneous Material Support:** Calculates heat transfer across boundaries with vastly different thermal conductivities (e.g., FR-4 insulator vs. Copper traces).
* **Real PCB Topology Parsing:** Integrates `pcb-tools` to parse production Gerber RS-274X files (`.gtl`/`.gbr`), automatically mapping copper polygons into a spatial thermal conductivity matrix $K(x,y)$.
* **Dual Solving Engines:**
    * *Transient Solver:* Simulates real-time heat propagation using an explicit Euler scheme strictly optimized via NumPy vectorization. Supports dynamic power profiles (e.g., PWM or component shutdown).
    * *Steady-State Solver:* Instantly calculates thermal equilibrium ($t \to \infty$) using the Jacobi iteration method.
* **Virtual Thermal Probes:** Places virtual thermocouples on critical components to log and plot temperature evolution over time (heating/cooling curves).
* **GUI & Multithreading:** Built with `PyQt6` and embedded `matplotlib`. CPU-bound mathematical operations are offloaded to background threads to maintain a responsive interface.

## 🧮 Mathematical Foundation

The core engine numerically solves the divergent form of the 2D Heat Equation to accurately model heat flux across heterogeneous materials (copper vs. substrate):

$$\rho c_p \frac{\partial u}{\partial t} = \nabla \cdot (K(x,y) \nabla u) + Q(x,y) - \frac{2h(u - T_{amb})}{d}$$

* **Convective Boundary Conditions:** Implemented using the "Ghost Node" (fictitious node) method to enforce Newton's Law of Cooling at the PCB edges, maintaining $O(\Delta x^2)$ spatial accuracy without slow `for` loops.
* **Numerical Stability:** The explicit transient solver strictly adheres to the 2D Courant–Friedrichs–Lewy (CFL) stability condition: $dt \le \frac{dx^2}{4\alpha}$.

## 🛠️ Installation & Usage

### Prerequisites
* Python 3.9+
* Required libraries: `numpy`, `matplotlib`, `PyQt6`, `pcb-tools`, `pycairo`

# Install dependencies
pip install -r requirements.txt

# Run the application
python main.py
Quick Start
Load Components: Edit or load a CSV file containing SMT component coordinates and power dissipation (Watts).

Load Topology (Optional): Click "Load Top Copper (Gerber)" to import a real .gtl file for accurate thermal spreading.

Run Simulation: Select either Steady-State for instant equilibrium results or Transient to watch the thermal wave propagate.

📂 Project Architecture
solver.py: The high-performance mathematical engine (Transient & Jacobi solvers, Matrix operations).

data_loader.py: Handles CSV parsing and Gerber RS-274X rasterization into physical index maps.

gui.py: PyQt6 application interface and multithreading architecture.

visualization.py: Matplotlib canvas rendering and virtual probe plotting.

config.py: Thermodynamic constants and solver parameters.

👨‍💻 About The Author
Developed as an R&D project to bridge the gap between hardware engineering, thermodynamic physics, and high-performance computational Python.

## License

Distributed under the MIT License.


