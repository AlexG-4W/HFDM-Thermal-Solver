"""
BGA Test Case Generator for HFDM-Thermal-Solver.

Generates:
1. A components CSV with a realistic BGA package + surrounding passives.
2. A heatsinks CSV with a heatsink on the BGA.
3. A synthetic Gerber file (.gbr) with copper pads, traces, and a ground pour —
   usable by load_gerber_to_k_matrix() for thermal conductivity mapping.

Output files are written to the `tests/bga_test/` directory.
"""

import os
import csv
import math

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "tests", "bga_test")


def generate_components_csv(filepath):
    """
    Generates a components_power.csv for a BGA test scenario.

    Layout (100x100 mm board):
      - U1: BGA-256 FPGA at center (50, 50), 35x35 mm, 8.0 W
      - U2: DDR memory (75, 70), 12x8 mm, 1.5 W
      - U3: Power regulator (15, 50), 5x5 mm, 3.0 W
      - C1-C4: Decoupling caps around BGA, 0.01 W each
    """
    components = [
        {"Designator": "U1_BGA", "Center_X_mm": 50.0, "Center_Y_mm": 50.0,
         "Width_mm": 35.0, "Length_mm": 35.0, "Power_Watts": 8.0},
        {"Designator": "U2_DDR", "Center_X_mm": 75.0, "Center_Y_mm": 70.0,
         "Width_mm": 12.0, "Length_mm": 8.0, "Power_Watts": 1.5},
        {"Designator": "U3_VREG", "Center_X_mm": 15.0, "Center_Y_mm": 50.0,
         "Width_mm": 5.0, "Length_mm": 5.0, "Power_Watts": 3.0},
        {"Designator": "C1", "Center_X_mm": 35.0, "Center_Y_mm": 35.0,
         "Width_mm": 2.0, "Length_mm": 1.2, "Power_Watts": 0.01},
        {"Designator": "C2", "Center_X_mm": 65.0, "Center_Y_mm": 35.0,
         "Width_mm": 2.0, "Length_mm": 1.2, "Power_Watts": 0.01},
        {"Designator": "C3", "Center_X_mm": 35.0, "Center_Y_mm": 65.0,
         "Width_mm": 2.0, "Length_mm": 1.2, "Power_Watts": 0.01},
        {"Designator": "C4", "Center_X_mm": 65.0, "Center_Y_mm": 65.0,
         "Width_mm": 2.0, "Length_mm": 1.2, "Power_Watts": 0.01},
    ]

    fieldnames = ["Designator", "Center_X_mm", "Center_Y_mm",
                  "Width_mm", "Length_mm", "Power_Watts"]

    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(components)

    print(f"  Components CSV: {filepath} ({len(components)} components)")
    return components


def generate_heatsinks_csv(filepath):
    """
    Generates a heatsinks.csv with a heatsink over the BGA.
    """
    heatsinks = [
        {"Designator": "HS1_BGA", "Center_X_mm": 50.0, "Center_Y_mm": 50.0,
         "Width_mm": 40.0, "Length_mm": 40.0, "Convection_H": 120.0},
        {"Designator": "HS2_VREG", "Center_X_mm": 15.0, "Center_Y_mm": 50.0,
         "Width_mm": 10.0, "Length_mm": 10.0, "Convection_H": 80.0},
    ]

    fieldnames = ["Designator", "Center_X_mm", "Center_Y_mm",
                  "Width_mm", "Length_mm", "Convection_H"]

    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(heatsinks)

    print(f"  Heatsinks CSV:  {filepath} ({len(heatsinks)} heatsinks)")


def generate_gerber(filepath, board_w_mm=100.0, board_h_mm=100.0):
    lines = []
    lines.append("%FSLAX25Y25*%")
    lines.append("%MOMM*%")
    lines.append("%ADD10C,0.5000*%")
    lines.append("%ADD11C,2.0000*%")
    lines.append("%ADD12C,0.2000*%")
    
    def coord(x_mm, y_mm):
        x_int = int(round(x_mm * 100000))
        y_int = int(round(y_mm * 100000))
        return f"X{x_int}Y{y_int}"

    lines.append("G01*")

    lines.append("D10*") 
    bga_cx, bga_cy = 50.0, 50.0
    pitch = 1.0
    x_start = bga_cx - (15) * pitch / 2
    y_start = bga_cy - (15) * pitch / 2
    for row in range(16):
        for col in range(16):
            px = x_start + col * pitch
            py = y_start + row * pitch
            lines.append(f"{coord(px, py)}D03*")

    lines.append("D11*")
    lines.append(f"{coord(15.0, 50.0)}D02*")
    lines.append(f"{coord(32.5, 50.0)}D01*")

    lines.append("D12*")
    for i in range(4):
        trace_y = 66.0 + i * 2.0
        lines.append(f"{coord(75.0, trace_y)}D02*")
        lines.append(f"{coord(57.5, trace_y)}D01*")
        
    lines.append("M02*")

    with open(filepath, "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=== BGA Test Case Generator ===\n")

    comp_path = os.path.join(OUTPUT_DIR, "bga_components.csv")
    hs_path = os.path.join(OUTPUT_DIR, "bga_heatsinks.csv")
    gerber_path = os.path.join(OUTPUT_DIR, "bga_top_copper.gbr")

    generate_components_csv(comp_path)
    generate_heatsinks_csv(hs_path)
    generate_gerber(gerber_path)

    print(f"\nAll files written to: {OUTPUT_DIR}")
    print("Done.")


if __name__ == "__main__":
    main()