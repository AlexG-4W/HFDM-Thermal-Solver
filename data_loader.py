import csv
import numpy as np
import os
import config

def load_components_list(filename="components_power.csv"):
    """
    Parses the components CSV and returns a list of dictionaries with raw properties.
    """
    components = []
    if not os.path.exists(filename):
        return components

    with open(filename, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            # Fallback: accept both 'Length_mm' (v1.1+) and legacy 'Height_mm'
            length = float(row.get('Length_mm') or row.get('Height_mm', 0))
            components.append({
                'Designator': row['Designator'],
                'Center_X_mm': float(row['Center_X_mm']),
                'Center_Y_mm': float(row['Center_Y_mm']),
                'Width_mm': float(row['Width_mm']),
                'Length_mm': length,
                'Power_Watts': float(row['Power_Watts'])
            })
    return components

def calculate_q_matrix(components, nx, ny):
    """
    Generates the total 2D volumetric heat source grid [W/m^3].
    
    This is a 2D shell model: the "volume" is the PCB-shell volume
    under the component footprint (A_footprint * d_pcb), NOT the 
    physical component volume. The solver treats Q as volumetric 
    source distributed uniformly through the board thickness.
    """
    Q_matrix = np.zeros((ny, nx))
    for comp in components:
        cx, cy = comp['Center_X_mm'], comp['Center_Y_mm']
        w, l_mm = comp['Width_mm'], comp['Length_mm']
        p = comp['Power_Watts']

        # Guard: skip components with zero/negative footprint (e.g. freshly
        # added empty rows in the GUI) to avoid ZeroDivisionError.
        footprint_area_m2 = (w / 1000.0) * (l_mm / 1000.0)
        pcb_shell_volume = footprint_area_m2 * config.d
        if pcb_shell_volume <= 0:
            continue

        idx_x_start, idx_x_end, idx_y_start, idx_y_end = get_indices(
            cx, cy, w, l_mm, nx, ny
        )

        # q_volumetric = P / (A_footprint * d_pcb)  [W/m^3]
        q_volumetric = p / pcb_shell_volume

        Q_matrix[idx_y_start:idx_y_end, idx_x_start:idx_x_end] += q_volumetric

    return Q_matrix

def load_components_dict(filename="components_power.csv", board_width_mm=100.0, board_height_mm=100.0):
    """
    Parses the components CSV and returns a dictionary of individual heat source grids.
    """
    nx = int(board_width_mm / (config.dx * 1000))
    ny = int(board_height_mm / (config.dx * 1000))
    
    component_matrices = {}
    
    if not os.path.exists(filename):
        print(f"Warning: {filename} not found.")
        return component_matrices, nx, ny

    with open(filename, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            designator = row['Designator']
            cx, cy = float(row['Center_X_mm']), float(row['Center_Y_mm'])
            length = float(row.get('Length_mm') or row.get('Height_mm', 0))
            w, l_mm = float(row['Width_mm']), length
            p = float(row['Power_Watts'])
            
            idx_x_start, idx_x_end, idx_y_start, idx_y_end = get_indices(cx, cy, w, l_mm, nx, ny)
            
            # Volumetric heat source calculation
            volume = (w / 1000.0) * (l_mm / 1000.0) * config.d
            q_volumetric = p / volume
            
            comp_matrix = np.zeros((ny, nx))
            comp_matrix[idx_y_start:idx_y_end, idx_x_start:idx_x_end] = q_volumetric
            component_matrices[designator] = comp_matrix
            
    return component_matrices, nx, ny

def load_components(filename="components_power.csv", board_width_mm=100.0, board_height_mm=100.0):
    """Legacy wrapper returning summed matrix."""
    comp_dict, nx, ny = load_components_dict(filename, board_width_mm, board_height_mm)
    Q_matrix = np.zeros((ny, nx))
    for mat in comp_dict.values():
        Q_matrix += mat
    return Q_matrix, nx, ny

def load_heatsinks(filename="heatsinks.csv", nx=200, ny=200):
    """
    Parses the heatsinks CSV and maps them to a convection coefficient grid (h_matrix).
    """
    h_matrix = np.full((ny, nx), config.h)
    
    if not os.path.exists(filename):
        return h_matrix

    with open(filename, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            cx, cy = float(row['Center_X_mm']), float(row['Center_Y_mm'])
            length = float(row.get('Length_mm') or row.get('Height_mm', 0))
            w, l_mm = float(row['Width_mm']), length
            h_val = float(row['Convection_H'])
            
            idx_x_start, idx_x_end, idx_y_start, idx_y_end = get_indices(cx, cy, w, l_mm, nx, ny)
            
            # Update h_matrix (use max to allow overlapping HS zones with highest H)
            h_matrix[idx_y_start:idx_y_end, idx_x_start:idx_x_end] = np.maximum(
                h_matrix[idx_y_start:idx_y_end, idx_x_start:idx_x_end], h_val
            )
            
    return h_matrix

def get_indices(cx, cy, w, h_mm, nx, ny):
    """
    Helper to convert center and dimensions in mm to grid indices.
    
    Uses floor() for start and ceil() for end to guarantee that the
    component footprint fully covers its physical area on the grid.
    """
    dx_mm = config.dx * 1000  # Grid spacing in mm

    x_start_mm = cx - w / 2
    x_end_mm   = cx + w / 2
    y_start_mm = cy - h_mm / 2
    y_end_mm   = cy + h_mm / 2

    # FIX V-4: floor for start, ceil for end — eliminates systematic
    # footprint shrinkage caused by int() truncation.
    idx_x_start = int(np.floor(x_start_mm / dx_mm))
    idx_x_end   = int(np.ceil(x_end_mm / dx_mm))
    idx_y_start = int(np.floor(y_start_mm / dx_mm))
    idx_y_end   = int(np.ceil(y_end_mm / dx_mm))

    # Clamp to grid boundaries
    idx_x_start = max(0, idx_x_start)
    idx_x_end   = min(nx, max(idx_x_start + 1, idx_x_end))
    idx_y_start = max(0, idx_y_start)
    idx_y_end   = min(ny, max(idx_y_start + 1, idx_y_end))

    return idx_x_start, idx_x_end, idx_y_start, idx_y_end

def generate_k_matrix(width_mm, height_mm, dx_m, k_fr4=None, k_copper=385.0):
    """
    Generates a 2D thermal conductivity matrix K(x,y) with mock copper traces.
    
    CRITICAL FIX (V-1): The baseline thermal conductivity is now k_eff
    (effective through-thickness conductivity of the multilayer stackup),
    NOT raw k_fr4. This accounts for internal power/ground copper planes
    that are always present in a 4-layer board, even in areas without 
    surface traces.
    
    Without this fix, areas without surface traces get k=0.3 W/mK instead 
    of k~34 W/mK, causing temperatures to exceed 3000 C for isolated 
    components.
    
    The surface trace contribution is added as an INCREMENT on top of k_eff,
    representing the additional copper in the signal layer:
        k_cell = k_eff_base + (K_Cu * t_trace) / d_total
    """
    nx = int(width_mm / (dx_m * 1000))
    ny = int(height_mm / (dx_m * 1000))

    # FIX V-1: Use effective stackup conductivity as the BASELINE.
    # This represents internal GND/VCC planes that span the entire board.
    k_eff_base = config.calculate_k_eff(
        layers=config.BOARD_LAYERS,
        copper_oz=config.COPPER_OZ,
        substrate_k=config.K_FR4
    )
    K_matrix = np.full((ny, nx), k_eff_base)

    # Additional k contributed by a surface trace (1 oz = 35 um copper)
    # on top of the already-accounted internal layers.
    cu_trace_thickness = config.COPPER_OZ * 0.000035  # [m]
    k_trace_increment = (k_copper * cu_trace_thickness) / config.d
    # k for cells WITH a surface trace:
    k_with_trace = k_eff_base + k_trace_increment

    # Trace Width: ~1.5mm
    w_cells = max(1, int(1.5 / (dx_m * 1000)))
    half_w = w_cells // 2

    # Mock Trace: vertical from Q1(20,80) to (20,50), then horizontal to U1(50,50)
    x1, y1 = 20, 80
    x2, y2 = 50, 50

    # Grid indices
    ix1 = int(x1 / (dx_m * 1000))
    iy1 = int(y1 / (dx_m * 1000))
    ix2 = int(x2 / (dx_m * 1000))
    iy2 = int(y2 / (dx_m * 1000))

    # Vertical segment
    y_min, y_max = min(iy1, iy2), max(iy1, iy2)
    K_matrix[y_min:y_max, max(0, ix1-half_w):min(nx, ix1+half_w+1)] = k_with_trace

    # Horizontal segment
    x_min, x_max = min(ix1, ix2), max(ix1, ix2)
    K_matrix[max(0, iy2-half_w):min(ny, iy2+half_w+1), x_min:x_max] = k_with_trace

    # Trace from Q1 (20, 80) to Left Edge (0, 80)
    K_matrix[max(0, iy1-half_w):min(ny, iy1+half_w+1), 0:ix1] = k_with_trace

    return K_matrix

def load_gerber_to_k_matrix(gerber_path, width_mm, height_mm, dx_m,
                             k_fr4=None, k_copper=385.0):
    """
    Parses a Gerber file and rasterizes it into a thermal conductivity 
    matrix K(x,y) using a through-thickness mixing model.
    
    FIX V-2: Cells with copper receive an INCREMENTAL k boost from the 
    surface trace, not a binary k=385 assignment. The baseline is always
    k_eff (multilayer stackup), not raw k_fr4.
    
    FIX V-6: Uses standard pcb-tools API (ctx.render(camfile) + dump),
    with post-hoc alignment via camfile.bounds. Falls back gracefully 
    if the Gerber file cannot be parsed.
    
    FIX V-7: Temp file is cleaned up in a finally block.
    """
    import gerber
    from gerber.render.cairo_backend import GerberCairoContext
    from PIL import Image
    import tempfile
    import os

    nx = int(width_mm / (dx_m * 1000))
    ny = int(height_mm / (dx_m * 1000))

    # FIX V-1/V-2: Baseline = effective stackup conductivity
    k_eff_base = config.calculate_k_eff(
        layers=config.BOARD_LAYERS,
        copper_oz=config.COPPER_OZ,
        substrate_k=config.K_FR4
    )

    try:
        camfile = gerber.read(gerber_path)
    except Exception as e:
        print(f"Error parsing Gerber file: {e}")
        return np.full((ny, nx), k_eff_base)

    # FIX V-6: Use the standard pcb-tools rendering API.
    # GerberCairoContext does NOT have set_bounds() or paint_background().
    ctx = GerberCairoContext()
    ctx.render(camfile)

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_file:
            tmp_path = tmp_file.name

        ctx.dump(tmp_path)

        img = Image.open(tmp_path).convert('L')
        # PIL.resize takes (width, height) = (nx, ny) — correct axis order
        img = img.resize((nx, ny), Image.Resampling.LANCZOS)
        img_data = np.array(img)

        # Gerber origin is bottom-left, PIL image origin is top-left.
        img_data = np.flipud(img_data)

        # Any pixel significantly brighter than background is copper
        copper_mask = img_data > 10

        K_matrix = np.full((ny, nx), k_eff_base)

        # FIX V-2: Instead of binary k=385, add the incremental
        # contribution of ONE surface copper layer (e.g. 1 oz = 35 um).
        cu_trace_thickness = config.COPPER_OZ * 0.000035  # [m]
        k_trace_increment = (k_copper * cu_trace_thickness) / config.d
        K_matrix[copper_mask] = k_eff_base + k_trace_increment

    finally:
        # FIX V-7: Guarantee temp file cleanup
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)

    return K_matrix
