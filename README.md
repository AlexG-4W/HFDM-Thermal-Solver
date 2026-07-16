# HFDM Thermal Solver ![Version](https://img.shields.io/badge/version-v1.2-blue)

A 2D Finite Difference Method (FDM) thermal simulator for Printed Circuit Boards (PCBs). This tool allows engineers to evaluate heat dissipation, test component placement, and analyze thermal conductivity using realistic board parameters and Gerber top copper layers.

## Features

- **Gerber Top Copper Parsing:** Accurately calculates thermal conductivity (K-matrix) based on traces, pads, and polygons using the modern **PyGerber** library (supports X2/X3 formats, precise rasterization).
- **Independent Data Input:** Seamlessly load Heat Sources (components) and Convection Zones (heatsinks) independently via CSV.
- **Steady-State & Transient Analysis:** Run both equilibrium thermal simulations and time-dependent temperature evolution.
- **Interactive Virtual Probes:** Left-click anywhere on the generated heatmap to place up to 10 virtual temperature probes. Right-click to clear them.
- **High-Quality Export:** Save the final thermal simulation results (heatmap + probes) directly to a PNG image for reporting.
- **Clean Executable Builds:** The build pipeline is optimized and free of legacy dependencies.


<img width="2545" height="1391" alt="scr11" src="https://github.com/user-attachments/assets/8ba90b6a-2f2f-4a94-9f6f-5557e183f6be" />



<img width="2551" height="1387" alt="scr22" src="https://github.com/user-attachments/assets/a73b86eb-a55f-475b-80bb-3702203f05c4" />






## Installation & Dependencies

Requires Python 3.9+

Install the required dependencies via `pip`:

```bash
pip install numpy matplotlib PyQt6 Pillow pygerber
```
*(Note: As of v1.2, the deprecated `pcb-tools` and `pycairo` libraries have been entirely removed and are no longer required).*

## Usage Workflow

The UI is built with PyQt6 and follows a straightforward workflow:

1. **Load Components CSV:** Load the power-dissipating components (Heat sources).
2. **Load Heatsinks CSV:** Load the local cooling zones (Convective heat transfer).
3. **Load Top Copper (Gerber):** Select a `.gbr` file to map copper traces and pads into the thermal conductivity matrix.
4. **Run Simulation:** Adjust ambient temperature (`T_amb`) and time (`t_final`), then click **Run Steady-State** or **Run Transient**.
5. **Probe & Save:** 
   - *Left-click* on the heatmap to place virtual temperature probes.
   - *Right-click* to clear them.
   - Use the UI save functionality to export the visualization to a PNG file.

## CSV Format Guide

Both components and heatsinks are imported via CSV files. **Important:** As of v1.2, the terminology has been updated to use `Length_mm` instead of `Height_mm` to avoid confusion with the physical 3D Z-axis in a 2D simulation.

### Components CSV (`components_power.csv`)
Defines the heat-generating parts on the PCB.

| Designator | Center_X_mm | Center_Y_mm | Width_mm | Length_mm | Power_Watts |
|------------|-------------|-------------|----------|-----------|-------------|
| U1_BGA     | 50.0        | 50.0        | 35.0     | 35.0      | 8.0         |
| U2_DDR     | 75.0        | 70.0        | 12.0     | 8.0       | 1.5         |

### Heatsinks CSV (`heatsinks.csv`)
Defines areas with enhanced convective cooling.

| Designator | Center_X_mm | Center_Y_mm | Width_mm | Length_mm | Convection_H |
|------------|-------------|-------------|----------|-----------|--------------|
| HS1_BGA    | 50.0        | 50.0        | 40.0     | 40.0      | 120.0        |
