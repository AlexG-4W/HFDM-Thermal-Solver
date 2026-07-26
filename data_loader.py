import csv
import io
import logging
import numpy as np
import os
import config
from PIL import Image
from pygerber.gerberx3.api import (
    ColorScheme,
    Rasterized2DLayer,
    Rasterized2DLayerParams,
)

logger = logging.getLogger(__name__)

# Thickness of one ounce of copper foil [m].
_COPPER_OZ_THICKNESS_M = 35e-6

# --- CSV plumbing (FIX H-2) --------------------------------------------------
# Spreadsheets are the real input format here, and they emit a BOM, semicolons
# and cp1251 depending on locale. Every one of those used to surface as a bare
# KeyError('Designator') with no hint of the actual cause.

_CSV_ENCODINGS = ("utf-8-sig", "cp1251")
_CSV_DELIMITERS = ",;\t"

COMPONENT_COLUMNS = ("Designator", "Center_X_mm", "Center_Y_mm",
                     "Width_mm", "Power_Watts")
HEATSINK_COLUMNS = ("Center_X_mm", "Center_Y_mm", "Width_mm", "Convection_H")

# Length is accepted under its legacy v1.0 name, so it is checked separately
# from the required-column list.
_LENGTH_ALIASES = ("Length_mm", "Height_mm")

# Key holding the source line number, used to make parse errors locatable.
_LINE_KEY = "__line__"


def _decode_csv(filename):
    """Reads a CSV as text, trying encodings in order. Returns (text, encoding)."""
    last_error = None
    for encoding in _CSV_ENCODINGS:
        try:
            with open(filename, mode="r", encoding=encoding, newline="") as handle:
                return handle.read(), encoding
        except UnicodeDecodeError as exc:
            last_error = exc
    raise ValueError(
        f"'{os.path.basename(filename)}' is not readable as text "
        f"(tried {', '.join(_CSV_ENCODINGS)}): {last_error}"
    )


def read_csv_rows(filename, required_columns=()):
    """
    Reads a CSV into a list of dicts, tolerating what spreadsheet tools emit.

    FIX H-2. Four separate failure modes are handled here:

    * ``utf-8-sig`` strips the byte-order mark Excel writes for "CSV UTF-8".
      Without it ``csv.DictReader`` reads the first field name as
      ``'\\ufeffDesignator'`` and every lookup raises KeyError.
    * The delimiter is sniffed, so ``;`` files from RU/EU Excel parse as well
      as comma files.
    * A cp1251 retry covers legacy files from older Windows tools.
    * Headers are validated up front, so a missing or misspelled column names
      itself and lists what was actually found, instead of surfacing as
      ``KeyError('Designator')``.

    Raises
    ------
    FileNotFoundError, ValueError
        Both carry the file name and enough context to act on.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"CSV file not found: {filename}")

    text, encoding = _decode_csv(filename)
    if not text.strip():
        raise ValueError(f"'{os.path.basename(filename)}' is empty.")

    try:
        delimiter = csv.Sniffer().sniff(
            text[:4096], delimiters=_CSV_DELIMITERS
        ).delimiter
    except csv.Error:
        delimiter = ","

    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    fieldnames = [name.strip() for name in (reader.fieldnames or [])]

    missing = [name for name in required_columns if name not in fieldnames]
    if missing:
        raise ValueError(
            f"'{os.path.basename(filename)}' is missing required column(s): "
            f"{', '.join(missing)}. Columns found: "
            f"{', '.join(fieldnames) if fieldnames else '(none)'}. "
            f"Detected delimiter '{delimiter}', encoding {encoding}."
        )

    rows = []
    for line_number, raw_row in enumerate(reader, start=2):
        row = {
            (key.strip() if isinstance(key, str) else key):
                (value.strip() if isinstance(value, str) else value)
            for key, value in raw_row.items()
        }
        if not any(row.get(name) for name in fieldnames):
            continue  # blank trailing line
        row[_LINE_KEY] = line_number
        rows.append(row)

    logger.info(
        "Read %d row(s) from '%s' (encoding %s, delimiter '%s').",
        len(rows), os.path.basename(filename), encoding, delimiter,
    )
    return rows


def _to_float(row, column, filename, default=None):
    """Parses one numeric cell, reporting file and line on failure. FIX H-2."""
    raw = row.get(column)
    text = "" if raw is None else str(raw).strip()
    if not text:
        if default is not None:
            return default
        raise ValueError(
            f"'{os.path.basename(filename)}' line {row.get(_LINE_KEY, '?')}: "
            f"column '{column}' is empty."
        )
    try:
        return float(text)
    except ValueError:
        pass
    # EU spreadsheets write 3,5 for 3.5 whenever the delimiter is ';'.
    if text.count(",") == 1 and "." not in text:
        try:
            return float(text.replace(",", "."))
        except ValueError:
            pass
    raise ValueError(
        f"'{os.path.basename(filename)}' line {row.get(_LINE_KEY, '?')}: "
        f"column '{column}' is not a number: {text!r}"
    )


def _length_of(row, filename):
    """Reads the length column under either its current or legacy name."""
    for alias in _LENGTH_ALIASES:
        if str(row.get(alias) or "").strip():
            return _to_float(row, alias, filename)
    raise ValueError(
        f"'{os.path.basename(filename)}' line {row.get(_LINE_KEY, '?')}: "
        f"needs a length column (one of {', '.join(_LENGTH_ALIASES)})."
    )

def load_components_list(filename="components_power.csv"):
    """
    Parses the components CSV and returns a list of dictionaries.

    FIX H-2: reading and validation go through :func:`read_csv_rows`, so BOM,
    delimiter and encoding variants all parse, and a malformed file raises an
    error naming the file, the line and the column. A missing file now raises
    FileNotFoundError instead of returning an empty list, which previously made
    a mistyped path look like an empty board.
    """
    rows = read_csv_rows(filename, COMPONENT_COLUMNS)
    components = []
    for row in rows:
        designator = str(row.get('Designator') or "").strip()
        if not designator:
            designator = f"C{len(components)}"
        components.append({
            'Designator': designator,
            'Center_X_mm': _to_float(row, 'Center_X_mm', filename),
            'Center_Y_mm': _to_float(row, 'Center_Y_mm', filename),
            'Width_mm': _to_float(row, 'Width_mm', filename),
            'Length_mm': _length_of(row, filename),
            # An empty power cell is legitimate (a passive part); a missing
            # Power_Watts COLUMN is not, and read_csv_rows already rejected it.
            'Power_Watts': _to_float(row, 'Power_Watts', filename, default=0.0),
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
        if w <= 0 or l_mm <= 0:
            logger.warning(
                "Component '%s' has a zero-area footprint (%.3f x %.3f mm); "
                "its %.3f W are NOT in the simulation.",
                comp.get('Designator', '?'), w, l_mm, p,
            )
            continue

        (idx_x_start, idx_x_end,
         idx_y_start, idx_y_end, raw_cells) = _footprint(cx, cy, w, l_mm, nx, ny)
        cells = (idx_x_end - idx_x_start) * (idx_y_end - idx_y_start)

        # FIX M-1: normalise by the DISCRETISED footprint, not the nominal one.
        # get_indices expands the footprint outward to whole cells (floor/ceil),
        # so P / (nominal_area * d) spread over those cells injects more than P
        # whenever the part is not grid-aligned - measured +66.7% for a
        # 2.0 x 1.2 mm capacitor on a 0.5 mm grid. Dividing by the un-clipped
        # cell count makes the grid total exactly P for an on-board part, and a
        # part hanging over the edge deposits only its on-board share.
        q_volumetric = p / (raw_cells * (config.dx ** 2) * config.d)

        # FIX H-3: say out loud how much power actually reaches the grid.
        reaching = q_volumetric * cells * (config.dx ** 2) * config.d
        if not _check_on_board(
            comp.get('Designator', '?'), cx, cy, nx, ny, cells, raw_cells,
            kind="Component",
            detail=f"{reaching:.3f} W of {p:.3f} W reach the grid.",
        ):
            continue

        Q_matrix[idx_y_start:idx_y_end, idx_x_start:idx_x_end] += q_volumetric

    return Q_matrix

def load_components_dict(filename="components_power.csv",
                         board_width_mm=None, board_height_mm=None):
    """
    Parses the components CSV and returns a dictionary of individual heat source grids.

    FIX L-2: the board size defaults to the configured one instead of a
    hardcoded 100 x 100 mm.
    """
    nx, ny = config.grid_size(board_width_mm, board_height_mm)


    component_matrices = {}

    # FIX L-4: print() raises AttributeError in a --windowed PyInstaller build,
    # where sys.stdout is None.
    if not os.path.exists(filename):
        logger.warning("Components file '%s' not found.", filename)
        return component_matrices, nx, ny

    for comp in load_components_list(filename):
        w, l_mm = comp['Width_mm'], comp['Length_mm']
        volume = (w / 1000.0) * (l_mm / 1000.0) * config.d
        if volume <= 0:
            logger.warning(
                "Component '%s' has a zero-area footprint (%.3f x %.3f mm); "
                "skipped.", comp['Designator'], w, l_mm,
            )
            continue

        (idx_x_start, idx_x_end, idx_y_start, idx_y_end, raw_cells) = _footprint(
            comp['Center_X_mm'], comp['Center_Y_mm'], w, l_mm, nx, ny
        )
        cells = (idx_x_end - idx_x_start) * (idx_y_end - idx_y_start)
        if not _check_on_board(
            comp['Designator'], comp['Center_X_mm'], comp['Center_Y_mm'],
            nx, ny, cells, raw_cells, kind="Component",
        ):
            continue

        comp_matrix = np.zeros((ny, nx))
        comp_matrix[idx_y_start:idx_y_end, idx_x_start:idx_x_end] = (
            comp['Power_Watts'] / volume
        )
        component_matrices[comp['Designator']] = comp_matrix

    return component_matrices, nx, ny

def load_components(filename="components_power.csv",
                    board_width_mm=None, board_height_mm=None):
    """Legacy wrapper returning summed matrix."""
    comp_dict, nx, ny = load_components_dict(filename, board_width_mm, board_height_mm)
    Q_matrix = np.zeros((ny, nx))
    for mat in comp_dict.values():
        Q_matrix += mat
    return Q_matrix, nx, ny

def load_heatsinks(filename="heatsinks.csv", nx=None, ny=None, base_h=None):
    """
    Parses the heatsinks CSV and maps them to a convection coefficient grid (h_matrix).

    FIX L-2: nx/ny default to the configured board grid rather than a hardcoded
    200 x 200, which silently assumed a 100 mm board at dx = 0.5 mm.

    base_h is the convection coefficient over the bare board; heatsink zones
    raise it locally. Passing it in rather than reading config.h keeps the
    caller in charge of the environment, the same way T_amb works.
    """
    if nx is None or ny is None:
        default_nx, default_ny = config.grid_size()
        nx = default_nx if nx is None else nx
        ny = default_ny if ny is None else ny
    base_h = config.h if base_h is None else float(base_h)
    h_matrix = np.full((ny, nx), base_h)

    # A missing heatsink file is a normal state, not an error: the board simply
    # has uniform convection. Component files are stricter (see
    # load_components_list), because there the absence means no heat sources.
    if not os.path.exists(filename):
        logger.info(
            "Heatsink file '%s' not found; using uniform h = %.1f W/m2K.",
            filename, base_h,
        )
        return h_matrix

    for row in read_csv_rows(filename, HEATSINK_COLUMNS):
        name = str(row.get('Designator') or "").strip() or f"HS@line{row[_LINE_KEY]}"
        cx = _to_float(row, 'Center_X_mm', filename)
        cy = _to_float(row, 'Center_Y_mm', filename)
        w = _to_float(row, 'Width_mm', filename)
        l_mm = _length_of(row, filename)
        h_val = _to_float(row, 'Convection_H', filename, default=base_h)

        (idx_x_start, idx_x_end,
         idx_y_start, idx_y_end, raw_cells) = _footprint(cx, cy, w, l_mm, nx, ny)
        cells = (idx_x_end - idx_x_start) * (idx_y_end - idx_y_start)
        if not _check_on_board(
            name, cx, cy, nx, ny, cells, raw_cells, kind="Heatsink",
        ):
            continue

        # Update h_matrix (use max to allow overlapping HS zones with highest H)
        h_matrix[idx_y_start:idx_y_end, idx_x_start:idx_x_end] = np.maximum(
            h_matrix[idx_y_start:idx_y_end, idx_x_start:idx_x_end], h_val
        )

    return h_matrix

def _footprint(cx, cy, w, h_mm, nx, ny):
    """
    Maps a footprint in mm to grid indices, plus its un-clipped cell count.

    Returns
    -------
    (idx_x_start, idx_x_end, idx_y_start, idx_y_end, raw_cells)
        The four indices are clamped to the grid and describe a possibly EMPTY
        slice. ``raw_cells`` is how many cells the footprint would occupy on an
        unbounded grid, so callers can tell "fully on board" from "clipped"
        from "entirely outside".

    FIX H-3: the previous clamping only bounded the start index from below and
    then forced ``end = max(start + 1, end)``. A component at cx = +150 mm
    produced the inverted slice ``x[290:200]`` - empty, so its power silently
    vanished - while one at cx = -50 mm was dragged onto the single edge column
    and injected there. Both ends are now clamped independently, an off-board
    footprint yields a genuinely empty slice, and the caller reports it.
    """
    dx_mm = config.dx * 1000  # Grid spacing in mm

    # Coordinates are relative to the board origin, not to the CAD sheet
    # origin. A board exported at Y = -45..-16 mm is entirely at negative Y in
    # sheet coordinates; without this shift every component on it reads as
    # off-board and its power is dropped.
    ox = config.BOARD_ORIGIN_X_MM
    oy = config.BOARD_ORIGIN_Y_MM

    # FIX V-4: floor for start, ceil for end — eliminates systematic
    # footprint shrinkage caused by int() truncation.
    raw_x_start = int(np.floor((cx - ox - w / 2) / dx_mm))
    raw_x_end   = int(np.ceil((cx - ox + w / 2) / dx_mm))
    raw_y_start = int(np.floor((cy - oy - h_mm / 2) / dx_mm))
    raw_y_end   = int(np.ceil((cy - oy + h_mm / 2) / dx_mm))

    # A sub-cell footprint still occupies one cell.
    if raw_x_end <= raw_x_start:
        raw_x_end = raw_x_start + 1
    if raw_y_end <= raw_y_start:
        raw_y_end = raw_y_start + 1
    raw_cells = (raw_x_end - raw_x_start) * (raw_y_end - raw_y_start)

    idx_x_start = min(nx, max(0, raw_x_start))
    idx_x_end   = min(nx, max(0, raw_x_end))
    idx_y_start = min(ny, max(0, raw_y_start))
    idx_y_end   = min(ny, max(0, raw_y_end))

    return idx_x_start, idx_x_end, idx_y_start, idx_y_end, raw_cells


def get_indices(cx, cy, w, h_mm, nx, ny):
    """
    Converts a centre and size in mm to grid indices.

    Thin wrapper over :func:`_footprint` kept for callers that do not need the
    out-of-bounds information. The returned slice may be empty when the
    footprint lies off the board.
    """
    return _footprint(cx, cy, w, h_mm, nx, ny)[:4]


def _check_on_board(name, cx, cy, nx, ny, cells, raw_cells, kind, detail=""):
    """
    Reports a footprint that misses the board, wholly or partly. FIX H-3.

    Returns True when at least one cell is on the board.
    """
    x0, y0, x1, y1 = config.board_extent()
    frame = f"board X {x0:.2f}..{x1:.2f}, Y {y0:.2f}..{y1:.2f} mm"
    if cells <= 0:
        logger.warning(
            "%s '%s' centred at (%.2f, %.2f) mm lies entirely outside the "
            "%s and is NOT part of the simulation.%s",
            kind, name, cx, cy, frame, f" {detail}" if detail else "",
        )
        return False
    if cells < raw_cells:
        logger.warning(
            "%s '%s' centred at (%.2f, %.2f) mm extends past the edge of the "
            "%s; only %.0f%% of its footprint is on the board.%s",
            kind, name, cx, cy, frame,
            100.0 * cells / raw_cells, f" {detail}" if detail else "",
        )
    return True

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

def _bounding_box_mm(rendering_result):
    """
    Extracts the Gerber bounding box in millimetres.

    PyGerber returns Offset objects whose as_millimeters() normalises both
    %MOMM*% and %MOIN*% sources, so the caller always gets millimetres.
    """
    box = rendering_result.get_properties().gerber_bounding_box
    return (
        float(box.min_x.as_millimeters()),
        float(box.min_y.as_millimeters()),
        float(box.max_x.as_millimeters()),
        float(box.max_y.as_millimeters()),
    )


def _copper_coverage(image):
    """
    Returns the copper mask as a single-channel PIL image.

    The alpha channel of COPPER_ALPHA is a clean {0, 255} mask. Luminance is
    NOT a coverage ramp - convert("L") maps the copper colour to a constant
    (100 for the default scheme), so any threshold above that silently blanks
    the whole layer. Alpha is used whenever it exists.
    """
    if "A" in image.getbands():
        return image.getchannel("A")
    return image.convert("L")


def load_gerber_to_k_matrix(gerber_path, nx, ny, dx_m, dy_m,
                            layers=None, copper_oz=None, dpmm=40.0):
    """
    Rasterises a Gerber copper layer into a thermal conductivity field.

    Returns
    -------
    numpy.ndarray, shape (ny, nx), dtype float64
        Through-thickness averaged conductivity [W/mK]. Row 0 is y = 0, so the
        result matches imshow(..., origin='lower').

    Physical model (FIX C-1)
    ------------------------
    Every cell is a through-thickness composite of the whole stackup, never
    bulk copper. The baseline is the full multilayer effective conductivity -
    internal GND/VCC planes span the entire board, so bare cells are NOT raw
    FR-4. A surface trace only adds its own foil, averaged over board
    thickness:

        k_base  = calculate_k_eff(BOARD_LAYERS, COPPER_OZ, K_FR4)
        k_delta = K_CU * copper_oz * 35um / d
        k_cell  = k_base + coverage * k_delta

    For a 4-layer 1 oz board this spans 33.96 .. 42.38 W/mK. Writing K_CU
    (385) into copper cells over a calculate_k_eff() baseline (0.30, because
    the parameterless call defaults to layers=0) produced a 1283:1 contrast
    instead of 1.25:1 and drove peak temperatures to 525 C.

    Geometry (FIX C-2)
    ------------------
    PyGerber rasterises the bounding box of the drawn apertures, not the board
    outline. That box is placed on the grid at its true millimetre coordinates
    and scaled by its true size, so copper registers against the component
    coordinates from the CSV. Resizing the raster straight to (nx, ny) instead
    stretched a 61 x 30 mm copper box across the whole 100 x 100 mm board -
    X by 1.66, Y by 3.40 - destroying registration entirely.

    Resampling (FIX M-5)
    --------------------
    Downsampling ~40 px/mm to 2 px/mm is a 20:1 reduction. NEAREST point-
    samples and deletes whole features: on the BGA fixture it dropped three of
    the four 0.2 mm DDR traces completely and zeroed 47.9% of copper-bearing
    cells, while inflating total copper area by 18.6%. BOX area-averages, so
    every trace survives and the alpha ramp becomes a genuine per-cell copper
    fraction feeding `coverage` above.

    Raises
    ------
    RuntimeError
        If the layer cannot be rendered or carries no usable geometry.
        Falling back to a stretched or uniform field is deliberately not done:
        that is the failure mode this function exists to remove, and a silently
        plausible wrong field is worse than a refusal.
    """
    raster, bbox_mm = rasterise_gerber(gerber_path, dpmm=dpmm)
    return project_copper(raster, bbox_mm, nx, ny, dx_m, dy_m,
                          layers=layers, copper_oz=copper_oz,
                          source=gerber_path)


def rasterise_gerber(gerber_path, dpmm=40.0):
    """
    Parses and renders a Gerber once.

    Returns
    -------
    (PIL.Image, (x0, y0, x1, y1))
        A single-channel copper mask and its physical extent in millimetres.

    Split out from the projection so that moving the board frame costs a cheap
    resample instead of re-parsing the file: on a real 60 x 33 mm board the
    parse is seconds, the projection is milliseconds.
    """
    try:
        layer = Rasterized2DLayer(
            options=Rasterized2DLayerParams(
                source_path=gerber_path,
                colors=ColorScheme.COPPER_ALPHA,
                dpi=int(round(dpmm * 25.4)),
            )
        )
        rendering_result = layer.render()
        raster = _copper_coverage(rendering_result.get_image())
    except Exception as exc:
        raise RuntimeError(
            f"Could not rasterise Gerber '{gerber_path}': {exc}"
        ) from exc

    try:
        bbox_mm = _bounding_box_mm(rendering_result)
    except (AttributeError, TypeError, ValueError) as exc:
        raise RuntimeError(
            f"Gerber '{gerber_path}' rendered but exposed no usable bounding "
            f"box ({type(exc).__name__}: {exc}); cannot register copper "
            f"against board coordinates."
        ) from exc

    logger.info(
        "Gerber '%s': %d x %d px raster over (%.2f, %.2f)..(%.2f, %.2f) mm.",
        os.path.basename(gerber_path), raster.width, raster.height, *bbox_mm,
    )
    return raster, bbox_mm


def project_copper(raster, bbox_mm, nx, ny, dx_m, dy_m,
                   layers=None, copper_oz=None, source="gerber"):
    """
    Places an already-rendered copper mask onto the board grid.

    Cheap enough to re-run whenever the board frame moves.
    """
    layers = config.BOARD_LAYERS if layers is None else layers
    copper_oz = config.COPPER_OZ if copper_oz is None else copper_oz

    if nx <= 0 or ny <= 0:
        raise ValueError(f"Grid must be positive, got nx={nx}, ny={ny}")
    if dx_m <= 0 or dy_m <= 0:
        raise ValueError(f"Cell size must be positive, got dx={dx_m}, dy={dy_m}")

    # --- Physical endpoints of the conductivity range (FIX C-1) -------------
    k_base = config.calculate_k_eff(
        layers=layers, copper_oz=copper_oz, substrate_k=config.K_FR4
    )
    k_delta = (config.K_CU * copper_oz * _COPPER_OZ_THICKNESS_M) / config.d
    k_matrix = np.full((ny, nx), k_base, dtype=np.float64)

    dx_mm = dx_m * 1000.0
    dy_mm = dy_m * 1000.0
    x0_mm, y0_mm, x1_mm, y1_mm = bbox_mm

    span_x_mm = x1_mm - x0_mm
    span_y_mm = y1_mm - y0_mm
    if span_x_mm <= 0 or span_y_mm <= 0:
        logger.warning(
            "Gerber '%s' contains no copper (bounding box %.3f x %.3f mm); "
            "returning uniform k = %.3f W/mK.",
            source, span_x_mm, span_y_mm, k_base,
        )
        return k_matrix

    # Positions are relative to the board origin: a board exported at negative
    # Y would otherwise project entirely off the grid and come back blank.
    ox_mm = config.BOARD_ORIGIN_X_MM
    oy_mm = config.BOARD_ORIGIN_Y_MM
    rx0, ry0 = x0_mm - ox_mm, y0_mm - oy_mm
    rx1, ry1 = x1_mm - ox_mm, y1_mm - oy_mm

    # Cell (i, j) covers x in [j*dx, (j+1)*dx], y in [i*dy, (i+1)*dy].
    # Snap OUTWARD to whole cells, then pad the raster to fill that aligned box
    # so one output pixel maps to exactly one grid cell. Rounding the bounding
    # box to the nearest cell instead shifts copper by up to half a cell, which
    # on a 0.5 mm grid is enough to drop a 0.2 mm trace into the wrong row.
    col0 = int(np.floor(rx0 / dx_mm))
    col1 = int(np.ceil(rx1 / dx_mm))
    row0 = int(np.floor(ry0 / dy_mm))
    row1 = int(np.ceil(ry1 / dy_mm))
    patch_w = max(1, col1 - col0)
    patch_h = max(1, row1 - row0)

    bx0, by0, bx1, by1 = config.board_extent()
    if x0_mm < bx0 - dx_mm or y0_mm < by0 - dy_mm \
            or x1_mm > bx1 + dx_mm or y1_mm > by1 + dy_mm:
        logger.warning(
            "Gerber copper spans (%.2f, %.2f)..(%.2f, %.2f) mm but the board "
            "frame is (%.2f, %.2f)..(%.2f, %.2f) mm; copper outside it will be "
            "clipped. Use 'Fit to data' if the frame is wrong.",
            x0_mm, y0_mm, x1_mm, y1_mm, bx0, by0, bx1, by1,
        )

    # --- Area-average resample to the exact cell footprint (FIX M-5) --------
    # Pad the raster out to the cell-aligned box before resampling; the pad is
    # bare laminate, so it is zero coverage.
    px_per_mm_x = raster.width / span_x_mm
    px_per_mm_y = raster.height / span_y_mm
    aligned = Image.new(
        "L",
        (max(1, int(np.ceil(patch_w * dx_mm * px_per_mm_x))),
         max(1, int(np.ceil(patch_h * dy_mm * px_per_mm_y)))),
        0,
    )
    # Raster row 0 is max_y, so the top inset is measured from the box's top.
    aligned.paste(
        raster,
        (int(round((rx0 - col0 * dx_mm) * px_per_mm_x)),
         int(round((row1 * dy_mm - ry1) * px_per_mm_y))),
    )
    patch = np.asarray(
        aligned.resize((patch_w, patch_h), Image.Resampling.BOX),
        dtype=np.float64,
    ) / 255.0
    # PyGerber puts max_y in row 0; the grid puts the board origin in row 0.
    patch = np.flipud(patch)
    np.clip(patch, 0.0, 1.0, out=patch)

    # --- Project onto the board, clipping whatever falls outside ------------
    dst_r0, dst_c0 = max(0, row0), max(0, col0)
    dst_r1 = min(ny, row0 + patch_h)
    dst_c1 = min(nx, col0 + patch_w)
    if dst_r1 <= dst_r0 or dst_c1 <= dst_c0:
        logger.warning(
            "Gerber copper lies entirely outside the board frame "
            "(%.2f, %.2f)..(%.2f, %.2f) mm; returning uniform k = %.3f W/mK.",
            bx0, by0, bx1, by1, k_base,
        )
        return k_matrix

    src_r0, src_c0 = dst_r0 - row0, dst_c0 - col0
    coverage = patch[src_r0:src_r0 + (dst_r1 - dst_r0),
                     src_c0:src_c0 + (dst_c1 - dst_c0)]

    k_matrix[dst_r0:dst_r1, dst_c0:dst_c1] = k_base + coverage * k_delta

    logger.info(
        "Copper projected onto cells [%d:%d, %d:%d]; mean coverage %.2f%%; "
        "k in [%.3f, %.3f] W/mK.",
        dst_r0, dst_r1, dst_c0, dst_c1,
        100.0 * coverage.mean(), k_matrix.min(), k_matrix.max(),
    )
    return k_matrix
