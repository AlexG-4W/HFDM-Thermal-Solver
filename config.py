# config.py
# Constants and Material Presets for 2D Thermal Simulation of a PCB

import logging

logger = logging.getLogger(__name__)

#: Single place the version is written down; the window title, the release
#: folder and the PyInstaller output all derive from it.
VERSION = "2.0"

# --- Physics Constants ---
# Copper Properties
K_CU = 385.0
RHO_CU = 8960.0
CP_CU = 390.0

# FR-4 Properties
K_FR4 = 0.3
RHO_FR4 = 1900.0
CP_FR4 = 1200.0

# --- Stackup Definition ---
d = 0.0016              # Total PCB thickness [m]
CU_THICKNESS_TOTAL = 0.00014  # Total copper thickness (e.g., 4 layers of 1oz/35um) [m]
FR4_THICKNESS = d - CU_THICKNESS_TOTAL

# --- Effective Properties Calculation (Parallel Mixing Rule) ---
# Effective Thermal Conductivity [W/mK]
k = (K_CU * CU_THICKNESS_TOTAL + K_FR4 * FR4_THICKNESS) / d

# Effective Density [kg/m^3]
rho = (RHO_CU * CU_THICKNESS_TOTAL + RHO_FR4 * FR4_THICKNESS) / d

# Effective Specific Heat [J/kgK]
# Note: cp is mass-weighted
cp = (CP_CU * RHO_CU * CU_THICKNESS_TOTAL + CP_FR4 * RHO_FR4 * FR4_THICKNESS) / (rho * d)

# --- Material Presets (Legacy support) ---
MATERIALS = {
    "FR-4": {"k": K_FR4, "rho": RHO_FR4, "cp": CP_FR4},
    "Aluminum": {"k": 200.0, "rho": 2700.0, "cp": 900.0},
    "Ceramic": {"k": 30.0, "rho": 3900.0, "cp": 800.0}
}

# --- Environment ---
T_amb = 25.0    # Ambient temperature [C]
h = 5.0         # Convective heat transfer coefficient [W/(m^2*K)]

# Surface emissivity for radiative loss. 0 reproduces the old convection-only
# model. Solder mask and bare FR-4 both sit near 0.9, and at 100 C radiation
# removes more heat than natural convection does, so 0 is not a safe default
# for a real board - it is only the backwards-compatible one for the library.
EMISSIVITY = 0.90

# --- Board Geometry (FIX L-2) ---
# The board size used to be the literal 100.0 repeated in ten places across
# gui.py, so a different board meant hunting for every one of them and the
# renderer and the solver could silently disagree about the extent.
BOARD_WIDTH_MM = 100.0
BOARD_HEIGHT_MM = 100.0

# Where the board sits in the CAD coordinate system. Real exports do not put
# the board at the sheet origin: an EasyEDA/Altium board commonly lands at
# X 31..95, Y -45..-16 mm. Without an origin every one of its components reads
# as off-board and its copper is projected outside the grid, which looks
# exactly like "the file failed to load".
BOARD_ORIGIN_X_MM = 0.0
BOARD_ORIGIN_Y_MM = 0.0

# --- Simulation Parameters ---
dx = 0.0005     # Grid spacing (0.5 mm) [m]
t_final = 150.0 # Total simulation time [s]

# Derived Quantities
alpha = k / (rho * cp)  # Thermal diffusivity [m^2/s]

# Maximum stable time step based on 2D CFL condition, with safety factor of 0.9
dt = 0.9 * (dx**2) / (4 * alpha)

# Global settings for reporting
SELECTED_MATERIAL = "FR-4 (Effective Multilayer)"
BOARD_LAYERS = 4
COPPER_OZ = 1

# Virtual Probes: (Name, X_mm, Y_mm)
probes = [
    ("U1_Core", 50.0, 50.0),
    ("Q1_Core", 20.0, 80.0),
    ("Edge_Sensor", 5.0, 5.0)
]

def update_derived_properties():
    """Recalculates alpha and dt based on current physical properties."""
    global alpha, dt
    alpha = k / (rho * cp)
    dt = 0.9 * (dx**2) / (4 * alpha)

def select_material(material_name):
    """
    Switches the substrate material and refreshes the derived quantities.

    FIX M-8: this was `pass`. Callers - including the test suite - changed
    nothing and got no error, so a test that selected a material and then
    asserted on the result was measuring the default state under another name.
    """
    global SELECTED_MATERIAL, k, rho, cp
    props = MATERIALS.get(material_name)
    if props is None:
        raise ValueError(
            f"Unknown material '{material_name}'. "
            f"Known materials: {', '.join(sorted(MATERIALS))}."
        )
    SELECTED_MATERIAL = material_name
    rho = props["rho"]
    cp = props["cp"]
    k = calculate_k_eff(layers=BOARD_LAYERS, copper_oz=COPPER_OZ,
                        substrate_k=props["k"])
    update_derived_properties()

def select_stackup(layers, copper_oz=1):
    """
    Switches the copper stackup and refreshes the derived quantities.

    FIX M-8: also a `pass` before. A test asking for layers=0 as its
    copper-free baseline silently kept BOARD_LAYERS=4, so the baseline and the
    4-layer case were the same board - which is why the comparison asserted
    26.805646 < 26.805646.
    """
    global BOARD_LAYERS, COPPER_OZ, k
    if layers < 0:
        raise ValueError(f"layers must be >= 0, got {layers}")
    if copper_oz <= 0:
        raise ValueError(f"copper_oz must be > 0, got {copper_oz}")
    BOARD_LAYERS = int(layers)
    COPPER_OZ = copper_oz
    substrate_k = MATERIALS.get(SELECTED_MATERIAL, {}).get("k", K_FR4)
    k = calculate_k_eff(layers=BOARD_LAYERS, copper_oz=COPPER_OZ,
                        substrate_k=substrate_k)
    update_derived_properties()

def board_extent():
    """(x0, y0, x1, y1) of the board in CAD millimetres."""
    return (BOARD_ORIGIN_X_MM,
            BOARD_ORIGIN_Y_MM,
            BOARD_ORIGIN_X_MM + BOARD_WIDTH_MM,
            BOARD_ORIGIN_Y_MM + BOARD_HEIGHT_MM)


def set_board_frame(origin_x_mm, origin_y_mm, width_mm, height_mm):
    """Moves and resizes the board frame in one place."""
    global BOARD_ORIGIN_X_MM, BOARD_ORIGIN_Y_MM, BOARD_WIDTH_MM, BOARD_HEIGHT_MM
    if width_mm <= 0 or height_mm <= 0:
        raise ValueError(
            f"Board must have positive size, got {width_mm} x {height_mm} mm"
        )
    BOARD_ORIGIN_X_MM = float(origin_x_mm)
    BOARD_ORIGIN_Y_MM = float(origin_y_mm)
    BOARD_WIDTH_MM = float(width_mm)
    BOARD_HEIGHT_MM = float(height_mm)


def grid_size(width_mm=None, height_mm=None, dx_m=None):
    """
    Number of cells (nx, ny) covering the board at the current cell size.

    FIX L-2: single source of truth for the grid dimensions. It also reports
    when the cell size does not tile the board exactly - at dx = 0.3 mm a
    100 mm board becomes 333 cells spanning 99.9 mm, and every coordinate on
    the plot is then off by that much with nothing to say so.
    """
    width_mm = BOARD_WIDTH_MM if width_mm is None else width_mm
    height_mm = BOARD_HEIGHT_MM if height_mm is None else height_mm
    dx_mm = (dx if dx_m is None else dx_m) * 1000.0
    if dx_mm <= 0:
        raise ValueError(f"Cell size must be positive, got {dx_mm} mm")

    nx = max(1, int(round(width_mm / dx_mm)))
    ny = max(1, int(round(height_mm / dx_mm)))
    for axis, cells, span in (("width", nx, width_mm), ("height", ny, height_mm)):
        covered = cells * dx_mm
        if abs(covered - span) > 1e-9:
            logger.warning(
                "Cell size %.4f mm does not tile the board %s of %.3f mm: "
                "%d cells cover %.3f mm (%+.3f mm).",
                dx_mm, axis, span, cells, covered, covered - span,
            )
    return nx, ny


def calculate_k_eff(layers=0, copper_oz=1, substrate_k=None):
    """
    Calculates the effective thermal conductivity of a multi-layer board using 
    the parallel mixing rule.
    """
    if substrate_k is None:
        substrate_k = K_FR4
        
    total_d = 0.0016
    cu_thickness_per_oz = 0.000035
    total_cu_d = layers * copper_oz * cu_thickness_per_oz
    
    if total_cu_d >= total_d:
        return K_CU
        
    fr4_d = total_d - total_cu_d
    k_eff = (K_CU * total_cu_d + substrate_k * fr4_d) / total_d
    return k_eff
