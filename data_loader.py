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
            components.append({
                'Designator': row['Designator'],
                'Center_X_mm': float(row['Center_X_mm']),
                'Center_Y_mm': float(row['Center_Y_mm']),
                'Width_mm': float(row['Width_mm']),
                'Height_mm': float(row['Height_mm']),
                'Power_Watts': float(row['Power_Watts'])
            })
    return components

def calculate_q_matrix(components, nx, ny):
    """
    Generates the total 2D heat source grid from a list of component dictionaries.
    """
    Q_matrix = np.zeros((ny, nx))
    for comp in components:
        cx, cy = comp['Center_X_mm'], comp['Center_Y_mm']
        w, h_mm = comp['Width_mm'], comp['Height_mm']
        p = comp['Power_Watts']
        
        idx_x_start, idx_x_end, idx_y_start, idx_y_end = get_indices(cx, cy, w, h_mm, nx, ny)
        
        # Volumetric heat source calculation
        volume = (w / 1000.0) * (h_mm / 1000.0) * config.d
        q_volumetric = p / volume
        
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
            w, h_mm = float(row['Width_mm']), float(row['Height_mm'])
            p = float(row['Power_Watts'])
            
            idx_x_start, idx_x_end, idx_y_start, idx_y_end = get_indices(cx, cy, w, h_mm, nx, ny)
            
            # Volumetric heat source calculation
            volume = (w / 1000.0) * (h_mm / 1000.0) * config.d
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
            w, h_mm = float(row['Width_mm']), float(row['Height_mm'])
            h_val = float(row['Convection_H'])
            
            idx_x_start, idx_x_end, idx_y_start, idx_y_end = get_indices(cx, cy, w, h_mm, nx, ny)
            
            # Update h_matrix (use max to allow overlapping HS zones with highest H)
            h_matrix[idx_y_start:idx_y_end, idx_x_start:idx_x_end] = np.maximum(
                h_matrix[idx_y_start:idx_y_end, idx_x_start:idx_x_end], h_val
            )
            
    return h_matrix

def get_indices(cx, cy, w, h_mm, nx, ny):
    """Helper to convert center and dimensions in mm to grid indices."""
    x_start_mm = cx - w/2
    x_end_mm = cx + w/2
    y_start_mm = cy - h_mm/2
    y_end_mm = cy + h_mm/2
    
    idx_x_start = int(x_start_mm / (config.dx * 1000))
    idx_x_end = int(x_end_mm / (config.dx * 1000))
    idx_y_start = int(y_start_mm / (config.dx * 1000))
    idx_y_end = int(y_end_mm / (config.dx * 1000))
    
    idx_x_start = max(0, idx_x_start)
    idx_x_end = min(nx, max(idx_x_start + 1, idx_x_end))
    idx_y_start = max(0, idx_y_start)
    idx_y_end = min(ny, max(idx_y_start + 1, idx_y_end))
    
    return idx_x_start, idx_x_end, idx_y_start, idx_y_end

def generate_k_matrix(width_mm, height_mm, dx_m, k_fr4=0.3, k_copper=385.0):
    """
    Generates a 2D thermal conductivity matrix K(x,y).
    Includes a mock copper trace connecting Q1 (20,80) and U1 (50,50).
    """
    nx = int(width_mm / (dx_m * 1000))
    ny = int(height_mm / (dx_m * 1000))
    K_matrix = np.full((ny, nx), k_fr4)
    
    # Trace Width: ~1.5mm (3 cells if dx=0.5mm)
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
    K_matrix[y_min:y_max, max(0, ix1-half_w):min(nx, ix1+half_w+1)] = k_copper
    
    # Horizontal segment
    x_min, x_max = min(ix1, ix2), max(ix1, ix2)
    K_matrix[max(0, iy2-half_w):min(ny, iy2+half_w+1), x_min:x_max] = k_copper
    
    # NEW: Trace to edge for cooling verification
    # Trace from Q1 (20, 80) to Left Edge (0, 80)
    K_matrix[max(0, iy1-half_w):min(ny, iy1+half_w+1), 0:ix1] = k_copper
    
    return K_matrix

def load_gerber_to_k_matrix(gerber_path, width_mm, height_mm, dx_m, k_fr4=0.3, k_copper=385.0):
    """
    Parses a Gerber file and rasterizes it into a thermal conductivity matrix K(x,y).
    """
    import gerber
    from gerber.render.cairo_backend import GerberCairoContext
    from PIL import Image
    import tempfile
    import os
    
    nx = int(width_mm / (dx_m * 1000))
    ny = int(height_mm / (dx_m * 1000))
    
    try:
        camfile = gerber.read(gerber_path)
    except Exception as e:
        print(f"Error parsing Gerber file: {e}")
        return np.full((ny, nx), k_fr4)
        
    ctx = GerberCairoContext()
    width_in = width_mm / 25.4
    height_in = height_mm / 25.4
    bounds_in = ((0.0, width_in), (0.0, height_in))
    
    # Enforce exact board dimensions to align traces perfectly with grid
    ctx.set_bounds(bounds_in)
    ctx.paint_background()
    ctx.new_render_layer()
    
    # Manually render primitives to bypass auto-bounding box
    for primitive in camfile.primitives:
        ctx.render(primitive)
        
    ctx.flatten()
    
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_file:
        tmp_path = tmp_file.name
        
    ctx.dump(tmp_path)
    
    img = Image.open(tmp_path).convert('L')
    # Resize the rendered image exactly to the simulation grid
    img = img.resize((nx, ny), Image.Resampling.LANCZOS)
    img_data = np.array(img)
    
    # Gerber origin is bottom-left, PIL image origin is top-left.
    # Our matplotlib plotting uses origin='lower' where index [0] corresponds to the bottom.
    img_data = np.flipud(img_data)
    
    # Any pixel significantly brighter than background (0) is considered copper
    copper_mask = img_data > 10
    
    K_matrix = np.full((ny, nx), k_fr4)
    K_matrix[copper_mask] = k_copper
    
    os.remove(tmp_path)
    
    return K_matrix
