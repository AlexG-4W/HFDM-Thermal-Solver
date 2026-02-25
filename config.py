# config.py
# Constants and Material Presets for 2D Thermal Simulation of a PCB

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
    """Legacy helper."""
    pass

def select_stackup(layers, copper_oz=1):
    """Legacy helper."""
    pass

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
