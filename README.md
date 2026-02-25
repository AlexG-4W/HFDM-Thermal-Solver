# HFDM Thermal Solver

A high-performance 2D Finite Difference Method (FDM) Thermal Solver designed to simulate transient and steady-state temperature distributions on Printed Circuit Boards (PCBs). 

The tool accounts for volumetric heat generation from SMT components, surface convective cooling from top and bottom faces, and lateral heat spreading through internal copper planes. It features a responsive PyQt6 desktop GUI, real-time interactive component editing, and native Gerber (RS-274X) file parsing for accurate physical modeling of copper topology.

## Features

- **Physics-Based Thermal Modeling**: Solves the 2D Heat Equation with a spatially variable thermal conductivity matrix $K(x,y)$.
- **Gerber Parsing Integration**: Load top copper (`.gbr`/`.gtl`) files to rasterize and map true trace geometries directly onto the conductivity matrix.
- **Interactive PyQt6 GUI**: A professional, multithreaded desktop application that stays responsive during intensive computations.
- **Multilayer Support**: Automatically calculates the effective thermal conductivity ($k_{eff}$) of 2-layer and 4-layer boards when specific trace data is unavailable.
- **Transient & Steady-State Modes**: Simulate temperature evolution over time or jump straight to the thermal equilibrium.
- **Interactive Component Table**: Edit component power dissipation, dimensions, and locations on the fly. 

## Installation

Ensure you have Python 3.8+ installed. Install the dependencies via `pip`:

```bash
pip install -r requirements.txt
```

*Note on Windows*: The `pycairo` library is required for rendering Gerber files.

## Usage

Launch the solver directly via Python:

```bash
python main.py
```

Alternatively, a standalone executable (`HFDM_Solver.exe`) is provided for Windows users, requiring no Python environment to run.

### GUI Workflow
1. **Load Data**: Start by loading a component definition CSV or manually adding components in the 'Components' tab.
2. **Load Top Copper (Optional)**: Click "Load Top Copper (Gerber)" to parse a physical PCB layout. Click "View Topology" to inspect the generated conductivity map.
3. **Configure Parameters**: Set the Ambient Temperature ($T_{amb}$) and the total time for transient simulations ($t_{final}$).
4. **Execute**: Run either the "Steady-State" or "Transient" simulation.
5. **Analyze**: View the live-updating heatmap indicating the temperature distribution across the board.

## File Structure

- `main.py`: Entry point for the application.
- `gui.py`: Main desktop application logic built with PyQt6.
- `solver.py`: Core computational FDM engine (Transient & Steady-State).
- `data_loader.py`: Handles CSV parsing and Gerber rasterization to build heat source ($Q$) and conductivity ($K$) matrices.
- `visualization.py`: Matplotlib plotting logic for heatmaps.
- `config.py`: Physical constants, stackup definitions, and global simulation parameters.
- `tests/`: Unit tests and physical validation scripts.

## License

Distributed under the MIT License.
