import sys
import logging
import os
import numpy as np
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QPushButton, QFileDialog, QGroupBox,
                             QLabel, QDoubleSpinBox, QTextEdit, QFrame, QSplitter,
                             QTableWidget, QTableWidgetItem, QTabWidget,
                             QRadioButton, QCheckBox, QButtonGroup, QSizePolicy)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QObject
from PyQt6.QtGui import QBrush, QColor

from board_view import BoardView
from data_loader import (load_heatsinks, load_components_list,
                         calculate_q_matrix, rasterise_gerber, project_copper)
from solver import PCBSolver
import config

# --- Logging Handler for GUI ---
class _LogBridge(QObject):
    """Carries log text from any thread into the GUI thread."""
    message = pyqtSignal(str)


class QTextEditHandler(logging.Handler):
    """
    Routes log records into the console pane.

    Records arrive from worker threads as well as the GUI thread - the solver
    logs its own convergence and interruption messages from inside run(). A
    QTextEdit may only be touched from the thread that owns it, so the text is
    handed over through a queued signal instead of being appended in place.
    Appending directly killed the process with STATUS_STACK_BUFFER_OVERRUN
    (0xC0000409) the first time the solver logged from the worker thread.
    """

    def __init__(self, text_edit):
        super().__init__()
        # Created on the GUI thread, so emitting from elsewhere queues the call.
        self._bridge = _LogBridge()
        self._bridge.message.connect(text_edit.append)

    def emit(self, record):
        try:
            self._bridge.message.emit(self.format(record))
        except RuntimeError:
            # The console widget has been destroyed; drop the record.
            pass

# --- Worker for Background Calculations ---
class SolverWorker(QObject):
    finished = pyqtSignal(object)  # Emits final u matrix
    progress = pyqtSignal(object, float, float) # Emits (u, current_time, max_temp)
    log = pyqtSignal(str)

    def __init__(self, nx, ny, Q, material, h_matrix, layers, copper_oz,
                 mode='steady', t_final=config.t_final, K_matrix=None,
                 T_amb=None, h_ambient=None, emissivity=0.0):
        super().__init__()
        self.nx = nx
        self.ny = ny
        self.Q = Q
        self.material = material
        self.h_matrix = h_matrix
        self.layers = layers
        self.copper_oz = copper_oz
        self.mode = mode
        self.t_final = t_final
        self.K_matrix = K_matrix
        # FIX L-7: the environment travels with the job instead of being read
        # off config from the worker thread while the GUI thread writes it.
        self.T_amb = T_amb
        self.h_ambient = h_ambient
        self.emissivity = emissivity
        self._is_running = True

    def stop(self):
        self._is_running = False

    def run(self):
        try:
            solver = PCBSolver(self.nx, self.ny, self.Q,
                               material_name=self.material,
                               h_matrix=self.h_matrix,
                               layers=self.layers,
                               copper_oz=self.copper_oz,
                               K_matrix=self.K_matrix,
                               T_amb=self.T_amb,
                               h_ambient=self.h_ambient,
                               emissivity=self.emissivity)
            
            self.log.emit(f"Starting {self.mode} simulation...")
            self.log.emit(f"Average k: {np.mean(solver.K):.4f} W/mK")
            if solver.emissivity > 0:
                self.log.emit(
                    f"Radiation on (eps = {solver.emissivity:.2f}); "
                    f"convection h = {np.min(solver.h_conv):.1f}"
                    f"..{np.max(solver.h_conv):.1f} W/m2K"
                )
            else:
                self.log.emit(
                    f"Radiation OFF; convection only, "
                    f"h = {np.min(solver.h_conv):.1f}..{np.max(solver.h_conv):.1f} W/m2K"
                )

            if self.mode == 'steady':
                # FIX C-5: hand the solver a cancellation probe so Stop and
                # window-close interrupt the sweep loop instead of waiting out
                # up to 50000 iterations.
                u_final, iterations, converged = solver.solve_steady_state(
                    should_continue=lambda: self._is_running
                )
                # FIX C-3: only claim convergence when it actually happened.
                if not self._is_running:
                    self.log.emit(
                        f"Steady-state STOPPED by user after {iterations} "
                        f"iterations - the result is not a solution."
                    )
                elif converged:
                    self.log.emit(f"Steady-state converged in {iterations} iterations.")
                else:
                    self.log.emit(
                        f"WARNING: steady-state did NOT converge in {iterations} "
                        f"iterations. Temperatures below are a lower bound, not "
                        f"a solution."
                    )
                self.finished.emit(u_final)

            elif self.mode == 'transient':
                total_steps = int(self.t_final / solver.dt)
                display_interval = 100
                
                for step in range(total_steps):
                    if not self._is_running:
                        break
                    
                    solver.step()
                    
                    if step % display_interval == 0:
                        t = step * solver.dt
                        max_t = np.max(solver.u)
                        self.progress.emit(np.copy(solver.u), t, max_t)
                
                self.log.emit("Transient simulation complete.")
                self.finished.emit(solver.u)

        except Exception as e:
            self.log.emit(f"Error: {str(e)}")
            self.finished.emit(None)

# --- Worker for Gerber rasterisation ---
class GerberWorker(QObject):
    """
    Runs load_gerber_to_k_matrix off the GUI thread.

    FIX H-4: the solver was already threaded but Gerber parsing was not. It ran
    inline in the click handler, preceded by QApplication.processEvents() - an
    antipattern that does not make the work asynchronous. It pumps the queue
    once BEFORE the blocking call, and allows re-entrancy: the user can trigger
    a second load from inside the nested event loop while the first is still
    rasterising.
    """
    finished = pyqtSignal(object)   # (raster, bbox_mm), or None on failure
    error = pyqtSignal(str)
    log = pyqtSignal(str)

    def __init__(self, path, dpmm=40.0):
        super().__init__()
        self.path = path
        self.dpmm = dpmm

    def run(self):
        try:
            self.log.emit(f"Rasterising {os.path.basename(self.path)}...")
            # Only the expensive half runs here. Placing the result on the
            # grid depends on the board frame, which the window may still
            # have to fit around this very file.
            self.finished.emit(rasterise_gerber(self.path, dpmm=self.dpmm))
        except Exception as exc:
            self.error.emit(f"Gerber Load Error: {exc}")
            self.finished.emit(None)


# --- Main Window ---
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"HFDM Thermal Solver v{config.VERSION}")
        self.resize(1200, 800)

        self.comp_dict = None
        self.components_list = []
        # Kept resident so switching layers never has to recompute it, and so
        # the renderer stops depending on which action last happened to build it.
        self.Q = None
        # Grid follows the configured board size and cell size (FIX L-2)
        self.nx, self.ny = config.grid_size()
        self.H = None
        self.K_matrix = None
        self.solver_thread = None
        self.worker = None
        self.gerber_thread = None
        self.gerber_worker = None
        self._gerber_path = None
        # The parsed raster is kept so moving the board frame is a cheap
        # reprojection instead of a fresh parse of the Gerber.
        self.gerber_raster = None
        self.gerber_bbox = None
        self.heatsink_path = None    # FIX M-2: remembered so H can be rebuilt

        # Virtual Probes state (up to 10)
        self.u_final = None          # Last completed temperature matrix
        # What is on screen lives in BoardView; probes read it from there so
        # the reading and the picture cannot drift apart (FIX M-4).
        self._has_result = False
        self._sim_phase = 'idle'     # 'idle' | 'running' | 'done'
        self._sim_mode = 'steady'
        self._sim_t = 0.0
        self._sim_tmax = None
        self.probes = []             # List of (x_mm, y_mm) probe positions
        self.probe_artists = []      # Matplotlib artists (marker + text) for cleanup
        self.MAX_PROBES = 10

        self.init_ui()
        self.setup_logging()
        self._sync_board_fields()
        self._sync_layer_controls()

    def _publish_probe_artists(self):
        """Hands the probe artists to the view so the Probes layer can hide them."""
        self.view.set_probe_artists(
            [artist for pair in self.probe_artists for artist in pair]
        )
        self._sync_layer_controls()

    # --- Layer controls -----------------------------------------------------
    # Base fields are mutually exclusive because there is one full-size image
    # and one colorbar; overlays are independent. The split is a consequence of
    # how Matplotlib composes a figure, not a stylistic choice.
    BASE_FIELD_CHOICES = (
        ('temperature', 'Temperature'),
        ('conductivity', 'Copper k'),
        ('power', 'Power density'),
        ('none', 'None'),
    )
    OVERLAY_CHOICES = (
        ('components', 'Component outlines'),
        ('copper', 'Copper mask'),
        ('convection', 'Heatsink zones'),
        ('probes', 'Probes'),
    )

    def _build_view_controls(self):
        """Layer panel, placed under the canvas rather than in the left rail.

        Switching layers is something you do while looking at the plot, so the
        control belongs in the plot's field of view; the left panel stays for
        loading data and setting parameters.
        """
        group = QGroupBox("View")
        # Two tight rows, never more. Without an explicit Fixed height policy
        # the group is Preferred like the canvas above it, so a maximised
        # window splits the spare vertical space between them and the strip
        # grows to half the pane with the rows drifting apart inside it.
        group.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        outer = QVBoxLayout()
        outer.setContentsMargins(8, 4, 8, 6)
        outer.setSpacing(4)

        base_row = QHBoxLayout()
        base_row.setSpacing(10)
        base_row.addWidget(QLabel("Base:"))
        self.base_group = QButtonGroup(self)
        self.radio_base = {}
        for key, text in self.BASE_FIELD_CHOICES:
            radio = QRadioButton(text)
            self.base_group.addButton(radio)
            self.radio_base[key] = radio
            radio.toggled.connect(
                lambda checked, k=key: self.on_base_field_changed(k) if checked else None
            )
            base_row.addWidget(radio)
        self.radio_base['temperature'].setChecked(True)
        base_row.addStretch()
        outer.addLayout(base_row)

        show_row = QHBoxLayout()
        show_row.setSpacing(10)
        show_row.addWidget(QLabel("Show:"))
        self.chk_layers = {}
        for key, text in self.OVERLAY_CHOICES:
            check = QCheckBox(text)
            check.setChecked(self.view.is_layer_visible(key))
            check.toggled.connect(
                lambda on, k=key: self.view.set_layer_visible(k, on)
            )
            self.chk_layers[key] = check
            show_row.addWidget(check)

        self.chk_lock_scale = QCheckBox("Lock colour scale")
        self.chk_lock_scale.setToolTip(
            "Hold the colour scale steady during a transient, so the same "
            "colour means the same temperature on every frame."
        )
        self.chk_lock_scale.toggled.connect(self.on_lock_scale_toggled)
        show_row.addWidget(self.chk_lock_scale)
        show_row.addStretch()
        outer.addLayout(show_row)

        group.setLayout(outer)
        return group

    # --- Board frame --------------------------------------------------------

    BOARD_FIT_MARGIN_MM = 1.0

    def data_extent(self):
        """
        Bounding box of everything loaded, in CAD millimetres, or None.

        Components and copper are quoted in the same CAD frame, so their union
        is the board as far as this program is concerned.
        """
        boxes = []
        for comp in self.components_list:
            w = comp.get('Width_mm') or 0.0
            length = comp.get('Length_mm') or 0.0
            if w <= 0 or length <= 0:
                continue
            boxes.append((comp['Center_X_mm'] - w / 2.0,
                          comp['Center_Y_mm'] - length / 2.0,
                          comp['Center_X_mm'] + w / 2.0,
                          comp['Center_Y_mm'] + length / 2.0))
        if self.gerber_bbox is not None:
            boxes.append(tuple(self.gerber_bbox))
        if not boxes:
            return None
        return (min(b[0] for b in boxes), min(b[1] for b in boxes),
                max(b[2] for b in boxes), max(b[3] for b in boxes))

    def _data_fits_frame(self):
        extent = self.data_extent()
        if extent is None:
            return True
        bx0, by0, bx1, by1 = config.board_extent()
        return (extent[0] >= bx0 and extent[1] >= by0
                and extent[2] <= bx1 and extent[3] <= by1)

    def fit_board_to_data(self, quiet=False):
        """Snaps the board frame around the loaded data, aligned to the grid."""
        extent = self.data_extent()
        if extent is None:
            if not quiet:
                self.log("Nothing loaded yet - nothing to fit the board to.")
            return False
        dx_mm = config.dx * 1000.0
        margin = self.BOARD_FIT_MARGIN_MM
        x0 = np.floor((extent[0] - margin) / dx_mm) * dx_mm
        y0 = np.floor((extent[1] - margin) / dx_mm) * dx_mm
        x1 = np.ceil((extent[2] + margin) / dx_mm) * dx_mm
        y1 = np.ceil((extent[3] + margin) / dx_mm) * dx_mm
        self.apply_board_frame(x0, y0, x1 - x0, y1 - y0)
        self.log(
            f"Board frame fitted to data: origin ({x0:.2f}, {y0:.2f}) mm, "
            f"{x1 - x0:.2f} x {y1 - y0:.2f} mm -> grid {self.nx} x {self.ny} "
            f"cells at dx = {dx_mm:.2f} mm."
        )
        return True

    def apply_board_frame(self, origin_x, origin_y, width, height):
        """
        Moves/resizes the board and rebuilds everything that depends on it.

        The grid, the convection field and the copper projection are all
        functions of the frame, so they are recomputed here rather than left
        for whoever notices first.
        """
        config.set_board_frame(origin_x, origin_y, width, height)
        self.nx, self.ny = config.grid_size()

        # The grid changed shape: the h field and the copper projection must
        # follow, and probe indices no longer mean anything.
        self.H = None
        self._ensure_h_matrix()
        self._project_gerber()
        self.clear_probes(redraw=False)
        self._has_result = False
        self.u_final = None
        if self.components_list:
            self.Q = calculate_q_matrix(self.components_list, self.nx, self.ny)
        self._sync_board_fields()
        self.show_layout()

    def _sync_board_fields(self):
        """Writes the current frame into the spin boxes without re-triggering."""
        values = {
            'origin_x': config.BOARD_ORIGIN_X_MM,
            'origin_y': config.BOARD_ORIGIN_Y_MM,
            'width': config.BOARD_WIDTH_MM,
            'height': config.BOARD_HEIGHT_MM,
        }
        for key, spin in self.board_spins.items():
            spin.blockSignals(True)
            spin.setValue(values[key])
            spin.blockSignals(False)
        self.lbl_grid.setText(f"Grid: {self.nx} x {self.ny} cells")

    def on_board_frame_edited(self, _value=None):
        self.apply_board_frame(
            self.board_spins['origin_x'].value(),
            self.board_spins['origin_y'].value(),
            self.board_spins['width'].value(),
            self.board_spins['height'].value(),
        )

    def _project_gerber(self):
        """Re-places the cached copper raster on the current grid (cheap)."""
        if self.gerber_raster is None or self.gerber_bbox is None:
            return
        try:
            self.K_matrix = project_copper(
                self.gerber_raster, self.gerber_bbox,
                self.nx, self.ny, config.dx, config.dx,
                source=os.path.basename(self._gerber_path or "gerber"),
            )
        except Exception as exc:
            self.logger.exception("Copper projection failed")
            self.log(f"Copper projection failed: {exc}")
            self.K_matrix = None
        self.view.set_copper(self.K_matrix)

    def on_convection_changed(self, _value=None):
        """Rebuilds the convection field when the base h changes."""
        self.H = None
        self._ensure_h_matrix()
        self._sync_layer_controls()
        self.log(
            f"Convection set to h = {self.spin_h.value():.2f} W/m2K"
            + (f" (heatsink zones up to {self.H.max():.0f})"
               if self.H is not None and self.H.max() > self.spin_h.value() else "")
        )

    def on_base_field_changed(self, kind):
        self.view.set_base_kind(kind)
        if kind != 'temperature':
            # The axes now carry W/mK or W/m3. A temperature probe on them
            # would be reading the wrong quantity and labelling it "°C"
            # (FIX M-4), so the probes go with the field they belong to.
            self.clear_probes(redraw=False)
        self.view.set_title(self._compose_title())

    def on_lock_scale_toggled(self, locked):
        self.view.set_clim_lock(locked, floor=self.spin_tamb.value())

    def _sync_layer_controls(self):
        """
        A layer switch is enabled exactly when its data exists.

        Disabled-with-a-reason also does the discovery work: a greyed-out
        "Copper mask (load a Gerber file)" tells the user what to do next.
        """
        available = {
            'temperature': self.view.base_array('temperature') is not None,
            'conductivity': self.K_matrix is not None,
            'power': self.Q is not None,
            'none': True,
        }
        reasons = {
            'conductivity': "Load a Gerber file to enable",
            'power': "Load or enter components to enable",
        }
        for key, radio in self.radio_base.items():
            radio.setEnabled(available[key])
            radio.setToolTip("" if available[key] else reasons.get(key, ""))
            if not available[key] and radio.isChecked():
                self.radio_base['temperature'].setChecked(True)

        overlay_available = {
            'components': bool(self.components_list),
            'copper': self.view.has_copper(),
            'convection': self.view.has_convection_zones(),
            'probes': bool(self.probes),
        }
        overlay_reasons = {
            'copper': "Load a Gerber file to enable",
            'convection': "Load a heatsinks CSV with a non-uniform h to enable",
            'probes': "Click the board after a run to place a probe",
        }
        for key, check in self.chk_layers.items():
            check.setEnabled(overlay_available[key])
            check.setToolTip(
                "" if overlay_available[key] else overlay_reasons.get(key, "")
            )

    # --- Plot title ---------------------------------------------------------

    def _compose_title(self):
        """
        Single composer for the plot title.

        Three different methods used to set it and they overwrote each other;
        now it is derived from state like everything else on the canvas.
        """
        kind = self.view.field_kind()
        if kind == 'conductivity':
            head = "Copper topology"
            if self.K_matrix is not None:
                head += (f" - k {self.K_matrix.min():.2f}"
                         f"..{self.K_matrix.max():.2f} W/mK")
        elif kind == 'power':
            head = "Power density"
        elif kind == 'none':
            head = "Board outline"
        elif self._sim_phase == 'running':
            head = (f"{self._sim_mode.capitalize()} running"
                    + (f" - t = {self._sim_t:.1f} s" if self._sim_mode == 'transient' else ""))
            if self._sim_tmax is not None:
                head += f" - Tmax {self._sim_tmax:.1f} °C"
        elif self._sim_phase == 'done':
            head = f"{self._sim_mode.capitalize()} result"
            if self._sim_tmax is not None:
                head += f" - Tmax {self._sim_tmax:.1f} °C"
        else:
            head = "Component layout"

        parts = [head]
        if self.components_list:
            parts.append(self._power_summary())
        return "  |  ".join(parts)

    def _power_summary(self):
        """How much of the nameplate power actually lands on the board (H-3)."""
        on_grid = (0.0 if self.Q is None
                   else float(self.Q.sum()) * (config.dx ** 2) * config.d)
        nameplate = sum(c['Power_Watts'] for c in self.components_list)
        text = f"{on_grid:.3f} W on the board"
        if abs(on_grid - nameplate) > 1e-9:
            text += f" of {nameplate:.3f} W declared"
        return text

    # BoardView owns the figure; these read-only views keep the probe and
    # export code readable without handing out a second owner of the canvas.
    @property
    def figure(self):
        return self.view.figure

    @property
    def canvas(self):
        return self.view.canvas

    @property
    def ax(self):
        return self.view.ax

    def init_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        # Splitter for adjustable panels
        splitter = QSplitter(Qt.Orientation.Horizontal)
        
        # --- Left Panel: Tabs ---
        left_panel = QFrame()
        left_panel.setMinimumWidth(400)
        left_layout = QVBoxLayout(left_panel)

        self.tabs = QTabWidget()
        
        # Tab 1: Setup
        setup_tab = QWidget()
        setup_layout = QVBoxLayout(setup_tab)

        # Data Group
        data_group = QGroupBox("Data Input")
        data_layout = QVBoxLayout()
        self.btn_load = QPushButton("Load Components CSV")
        self.btn_load.clicked.connect(self.load_data)
        data_layout.addWidget(self.btn_load)
        
        self.btn_save_csv = QPushButton("Save Edited CSV")
        self.btn_save_csv.clicked.connect(self.save_data_to_csv)
        self.btn_save_csv.setEnabled(False)
        data_layout.addWidget(self.btn_save_csv)

        self.lbl_status = QLabel("Status: No data loaded")
        data_layout.addWidget(self.lbl_status)
        
        self.btn_load_heatsinks = QPushButton("Load Heatsinks CSV")
        self.btn_load_heatsinks.clicked.connect(self.load_heatsinks_csv)
        data_layout.addWidget(self.btn_load_heatsinks)
        
        self.lbl_heatsinks_status = QLabel("Heatsinks: None")
        data_layout.addWidget(self.lbl_heatsinks_status)
        
        # Gerber Input
        self.btn_load_gerber = QPushButton("Load Top Copper (Gerber)")
        self.btn_load_gerber.clicked.connect(self.load_gerber)
        data_layout.addWidget(self.btn_load_gerber)
        
        self.lbl_gerber_status = QLabel("Gerber: None loaded")
        data_layout.addWidget(self.lbl_gerber_status)

        # "View Topology" is gone: it was a display MODE, and a mode you can
        # enter but not leave is exactly what made the canvas exclusive. Its
        # job is now the "Copper conductivity" base-field radio button.

        data_group.setLayout(data_layout)
        setup_layout.addWidget(data_group)

        # Parameters Group
        param_group = QGroupBox("Parameters")
        param_layout = QVBoxLayout()
        
        # Ambient Temp
        h_layout1 = QHBoxLayout()
        h_layout1.addWidget(QLabel("T_amb [°C]:"))
        self.spin_tamb = QDoubleSpinBox()
        self.spin_tamb.setRange(-50, 200)
        self.spin_tamb.setValue(config.T_amb)
        h_layout1.addWidget(self.spin_tamb)
        param_layout.addLayout(h_layout1)

        # t_final
        h_layout2 = QHBoxLayout()
        h_layout2.addWidget(QLabel("t_final [s]:"))
        self.spin_tfinal = QDoubleSpinBox()
        self.spin_tfinal.setRange(1, 3600)
        self.spin_tfinal.setValue(config.t_final)
        h_layout2.addWidget(self.spin_tfinal)
        param_layout.addLayout(h_layout2)

        # Base convection. Previously only reachable by editing config.py, so
        # the only way to model anything but still air was a dummy heatsink
        # covering the whole board.
        h_layout3 = QHBoxLayout()
        h_layout3.addWidget(QLabel("h [W/m²K]:"))
        self.spin_h = QDoubleSpinBox()
        self.spin_h.setRange(0.1, 10000.0)
        self.spin_h.setDecimals(2)
        self.spin_h.setValue(config.h)
        self.spin_h.setKeyboardTracking(False)
        self.spin_h.setToolTip(
            "Convection coefficient over the whole board. Still air 5-12, "
            "gentle airflow 20-30, forced air 40-80. Heatsink zones from the "
            "CSV override it locally."
        )
        self.spin_h.valueChanged.connect(self.on_convection_changed)
        h_layout3.addWidget(self.spin_h)
        param_layout.addLayout(h_layout3)

        h_layout4 = QHBoxLayout()
        h_layout4.addWidget(QLabel("Emissivity ε:"))
        self.spin_emissivity = QDoubleSpinBox()
        self.spin_emissivity.setRange(0.0, 1.0)
        self.spin_emissivity.setDecimals(2)
        self.spin_emissivity.setSingleStep(0.05)
        self.spin_emissivity.setValue(config.EMISSIVITY)
        self.spin_emissivity.setKeyboardTracking(False)
        self.spin_emissivity.setToolTip(
            "Surface emissivity for radiative loss. 0 disables radiation; "
            "solder mask and FR-4 are about 0.9. At 100 °C radiation carries "
            "more heat than natural convection, so leaving it at 0 makes a "
            "board look far hotter than it is."
        )
        h_layout4.addWidget(self.spin_emissivity)
        param_layout.addLayout(h_layout4)

        # --- Board frame -----------------------------------------------------
        # A CAD export does not place the board at the sheet origin, so the
        # frame needs an origin as well as a size. Auto-fitted on load, and
        # editable for the case where only part of a board is of interest.
        self.board_spins = {}
        for key, label, lo, hi in (
            ('origin_x', "Origin X [mm]:", -100000.0, 100000.0),
            ('origin_y', "Origin Y [mm]:", -100000.0, 100000.0),
            ('width', "Width [mm]:", 0.5, 100000.0),
            ('height', "Height [mm]:", 0.5, 100000.0),
        ):
            row = QHBoxLayout()
            row.addWidget(QLabel(label))
            spin = QDoubleSpinBox()
            spin.setRange(lo, hi)
            spin.setDecimals(3)
            spin.setSingleStep(1.0)
            # Without this, every keystroke re-grids the board mid-typing.
            spin.setKeyboardTracking(False)
            spin.valueChanged.connect(self.on_board_frame_edited)
            self.board_spins[key] = spin
            row.addWidget(spin)
            param_layout.addLayout(row)

        fit_row = QHBoxLayout()
        self.lbl_grid = QLabel("Grid: -")
        fit_row.addWidget(self.lbl_grid)
        self.btn_fit_board = QPushButton("Fit to data")
        self.btn_fit_board.setToolTip(
            "Set the frame from the loaded components and copper."
        )
        self.btn_fit_board.clicked.connect(lambda: self.fit_board_to_data())
        fit_row.addWidget(self.btn_fit_board)
        param_layout.addLayout(fit_row)

        param_group.setLayout(param_layout)
        setup_layout.addWidget(param_group)

        # Execution Group
        exec_group = QGroupBox("Execution")
        exec_layout = QVBoxLayout()
        self.btn_steady = QPushButton("Run Steady-State")
        self.btn_steady.clicked.connect(lambda: self.start_simulation('steady'))
        self.btn_transient = QPushButton("Run Transient")
        self.btn_transient.clicked.connect(lambda: self.start_simulation('transient'))
        self.btn_stop = QPushButton("Stop Simulation")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self.stop_simulation)
        
        exec_layout.addWidget(self.btn_steady)
        exec_layout.addWidget(self.btn_transient)
        exec_layout.addWidget(self.btn_stop)
        exec_group.setLayout(exec_layout)
        setup_layout.addWidget(exec_group)
        setup_layout.addStretch()

        self.tabs.addTab(setup_tab, "Setup")

        # Tab 2: Components
        comp_tab = QWidget()
        comp_layout = QVBoxLayout(comp_tab)
        
        self.table = QTableWidget()
        self.table.setColumnCount(6)
        self.table.setHorizontalHeaderLabels(["Designator", "Power [W]", "X [mm]", "Y [mm]", "W [mm]", "L [mm]"])
        self.table.cellChanged.connect(self.on_table_edit)
        comp_layout.addWidget(self.table)

        btn_comp_layout = QHBoxLayout()
        self.btn_add_row = QPushButton("Add Component")
        self.btn_add_row.clicked.connect(self.add_component_row)
        self.btn_del_row = QPushButton("Delete Selected")
        self.btn_del_row.clicked.connect(self.delete_component_row)
        btn_comp_layout.addWidget(self.btn_add_row)
        btn_comp_layout.addWidget(self.btn_del_row)
        comp_layout.addLayout(btn_comp_layout)
        
        self.tabs.addTab(comp_tab, "Components")
        left_layout.addWidget(self.tabs)

        # Log Console
        left_layout.addWidget(QLabel("Console Output:"))
        self.log_console = QTextEdit()
        self.log_console.setReadOnly(True)
        self.log_console.setStyleSheet("background-color: #1e1e1e; color: #d4d4d4; font-family: Consolas;")
        left_layout.addWidget(self.log_console)

        splitter.addWidget(left_panel)

        # --- Right Panel: Visualization ---
        right_panel = QFrame()
        right_layout = QVBoxLayout(right_panel)
        
        # The canvas and every artist on it belong to BoardView. MainWindow
        # owns data and threads, and pushes data in; it never draws.
        self.view = BoardView()
        # The canvas takes every spare pixel; the strips below keep their
        # natural height.
        right_layout.addWidget(self.view, 1)

        # Connect mouse click on canvas for interactive probes
        self.view.canvas.mpl_connect('button_press_event', self.on_canvas_click)

        right_layout.addWidget(self._build_view_controls(), 0)

        # Probe & Export toolbar under the canvas
        toolbar_layout = QHBoxLayout()
        self.btn_clear_probes = QPushButton("Clear Probes")
        # Wrapped: clicked emits checked=False, which would bind to the new
        # redraw parameter and clear the probes without repainting the canvas.
        self.btn_clear_probes.clicked.connect(lambda: self.clear_probes())
        toolbar_layout.addWidget(self.btn_clear_probes)

        self.btn_save_image = QPushButton("Save Result Image")
        self.btn_save_image.clicked.connect(self.save_result_image)
        toolbar_layout.addWidget(self.btn_save_image)
        right_layout.addLayout(toolbar_layout)

        splitter.addWidget(right_panel)
        
        main_layout.addWidget(splitter)

    # Loggers whose records are mirrored into the console pane. data_loader is
    # on the list because the H-3 out-of-bounds warnings are only useful if the
    # user actually sees them; attaching to the root logger instead would drag
    # in matplotlib/PIL/pygerber chatter.
    LOGGED_MODULES = ("HFDM_GUI", "data_loader", "solver")

    def setup_logging(self):
        self.handler = QTextEditHandler(self.log_console)
        self.handler.setFormatter(logging.Formatter('%(asctime)s - %(message)s', '%H:%M:%S'))
        for name in self.LOGGED_MODULES:
            module_logger = logging.getLogger(name)
            if self.handler not in module_logger.handlers:
                module_logger.addHandler(self.handler)
            module_logger.setLevel(logging.INFO)
        self.logger = logging.getLogger("HFDM_GUI")
        self.logger.info("HFDM GUI Initialized.")

    def log(self, message):
        self.logger.info(message)

    def load_data(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open CSV", "", "CSV Files (*.csv);;All Files (*)"
        )
        if not path:
            return

        try:
            self.components_list = load_components_list(path)
            # Ensure nx, ny are fresh based on current config
            self.nx, self.ny = config.grid_size()
            # FIX M-2: do NOT rebuild H here. Loading components must not touch
            # the convection layer.
            self._ensure_h_matrix()

            self.populate_table()
            if not self._data_fits_frame():
                self.log("Components fall outside the current board frame.")
                self.fit_board_to_data()
            
            self.lbl_status.setText(f"Loaded: {self.nx}x{self.ny} grid")
            self.log(f"Successfully loaded {len(self.components_list)} components and heatsinks.")
            self.btn_save_csv.setEnabled(True)
            
            # Show initial empty board or placeholder
            self.on_table_edit() # This triggers initial rendering
            
        except Exception as e:
            self.logger.exception("Components Load Error")
            self.log(f"Load Error: {str(e)}")

    def load_heatsinks_csv(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Heatsinks CSV", "", "CSV Files (*.csv);;All Files (*)"
        )
        if not path:
            return

        try:
            from data_loader import load_heatsinks, read_csv_rows, HEATSINK_COLUMNS

            # FIX H-2: count through the same validated reader instead of a
            # second bare open(), which ignored BOM, delimiter and encoding.
            count = len(read_csv_rows(path, HEATSINK_COLUMNS))
            self.H = load_heatsinks(path, nx=self.nx, ny=self.ny,
                                    base_h=self.spin_h.value())
            self.heatsink_path = path   # FIX M-2: survives a component reload
            self.view.set_convection(self.H)
            if self.view.has_convection_zones():
                self.chk_layers['convection'].setChecked(True)
            self._sync_layer_controls()

            self.lbl_heatsinks_status.setText(f"Heatsinks: {count} loaded")
            self.log(f"Successfully loaded {count} heatsinks.")

        except Exception as e:
            self.logger.exception("Heatsinks Load Error")
            self.log(f"Heatsinks Load Error: {str(e)}")

    # ----- Busy state (FIX H-5) -----

    # Only what MUTATES data the worker is reading. View controls are not here
    # on purpose - see the note in _set_busy.
    BUSY_WIDGETS = (
        "btn_load", "btn_save_csv", "btn_load_heatsinks", "btn_load_gerber",
        "btn_add_row", "btn_del_row", "btn_steady", "btn_transient",
        "btn_fit_board",
    )

    def _set_busy(self, busy, can_stop=False):
        """
        Locks the UI down for the duration of a background job.

        FIX H-5: only btn_steady, btn_transient and btn_stop used to change
        state, leaving eight controls and the component table live during a
        run. Editing a cell mid-transient rebuilt the whole figure, destroying
        the AxesImage that the worker's next progress signal then tried to
        update through a dangling reference.

        That lock is now narrower on purpose. The reason it had to cover the
        view was fig.clear(); with artists persistent, changing what is drawn
        no longer touches what the worker is driving, and both the click and
        the progress signal are handled on the GUI thread, so they serialise.
        Toggling the copper mask while watching a transient is a thing people
        want to do, so the layer panel, Clear Probes and Save Image stay live.
        The table and the loaders stay locked: they change data, not the view.

        can_stop is False for jobs that cannot be interrupted, so the Stop
        button never offers something the code cannot deliver: the Gerber
        rasteriser runs inside PyGerber and has no cancellation point.
        """
        for name in self.BUSY_WIDGETS:
            widget = getattr(self, name, None)
            if widget is not None:
                widget.setEnabled(not busy)
        self.table.setEnabled(not busy)
        self.spin_tamb.setEnabled(not busy)
        self.spin_tfinal.setEnabled(not busy)
        for spin in self.board_spins.values():
            spin.setEnabled(not busy)
        self.spin_h.setEnabled(not busy)
        self.spin_emissivity.setEnabled(not busy)
        self.btn_stop.setEnabled(busy and can_stop)
        if not busy:
            # Restore the CORRECT idle state, not merely "everything on":
            # Save CSV has its own precondition and must not come back enabled
            # just because a run ended.
            self.btn_save_csv.setEnabled(bool(self.components_list))

    # ----- Gerber loading (FIX H-4) -----

    def load_gerber(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Top Copper Gerber", "",
            "Gerber Files (*.gbr *.gtl);;All Files (*)"
        )
        if not path:
            return

        self._gerber_path = path
        self.lbl_gerber_status.setText("Gerber: loading...")
        self.log(f"Loading Gerber '{os.path.basename(path)}' in the background...")
        self._set_busy(True, can_stop=False)

        self.gerber_thread = QThread()
        self.gerber_worker = GerberWorker(path)
        self.gerber_worker.moveToThread(self.gerber_thread)

        self.gerber_thread.started.connect(self.gerber_worker.run)
        self.gerber_worker.finished.connect(self.on_gerber_loaded)
        self.gerber_worker.finished.connect(self.gerber_thread.quit)
        self.gerber_worker.finished.connect(self.gerber_worker.deleteLater)
        self.gerber_thread.finished.connect(self.gerber_thread.deleteLater)
        self.gerber_worker.log.connect(self.log)
        self.gerber_worker.error.connect(self.log)

        self.gerber_thread.start()

    def on_gerber_loaded(self, result):
        # References are deliberately not cleared here - see the note in
        # on_simulation_finished about destroying a worker mid-emission.
        if result is None:
            self.lbl_gerber_status.setText("Gerber: load failed")
        else:
            self.gerber_raster, self.gerber_bbox = result
            filename = os.path.basename(self._gerber_path or "")
            self.lbl_gerber_status.setText(f"Gerber: {filename}")
            x0, y0, x1, y1 = self.gerber_bbox
            self.log(
                f"Gerber '{filename}' covers ({x0:.2f}, {y0:.2f})..."
                f"({x1:.2f}, {y1:.2f}) mm."
            )
            # A CAD export rarely sits inside a default 0..100 mm frame; fit
            # the board around it rather than silently clipping it away.
            if self._data_fits_frame():
                self._project_gerber()
            else:
                self.log("Copper falls outside the current board frame.")
                self.fit_board_to_data()
            # Show it straight away: the whole point of loading a Gerber is to
            # check it against the component coordinates.
            self.chk_layers['copper'].setChecked(True)
            if self.K_matrix is not None:
                self.log(
                    f"Successfully loaded Gerber: {filename} "
                    f"(k from {self.K_matrix.min():.2f} to "
                    f"{self.K_matrix.max():.2f} W/mK)"
                )
        self._set_busy(False)
        self._sync_layer_controls()
        self.view.set_title(self._compose_title())

    def populate_table(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.components_list))
        for i, comp in enumerate(self.components_list):
            self.table.setItem(i, 0, QTableWidgetItem(str(comp['Designator'])))
            self.table.setItem(i, 1, QTableWidgetItem(str(comp['Power_Watts'])))
            self.table.setItem(i, 2, QTableWidgetItem(str(comp['Center_X_mm'])))
            self.table.setItem(i, 3, QTableWidgetItem(str(comp['Center_Y_mm'])))
            self.table.setItem(i, 4, QTableWidgetItem(str(comp['Width_mm'])))
            self.table.setItem(i, 5, QTableWidgetItem(str(comp['Length_mm'])))
        self.table.blockSignals(False)

    def _ensure_h_matrix(self):
        """
        Guarantees a convection field exists, without ever discarding one.

        FIX M-2: load_data() called load_heatsinks() with the DEFAULT filename
        on every component load, so a heatsink CSV the user had already
        selected was silently replaced by the uniform default. Components and
        heatsinks are independent input layers; neither may overwrite the
        other. The path is remembered so the field can be rebuilt if the grid
        ever changes, instead of being thrown away.
        """
        if self.H is not None and self.H.shape == (self.ny, self.nx):
            return
        if self.heatsink_path and os.path.exists(self.heatsink_path):
            self.H = load_heatsinks(self.heatsink_path, nx=self.nx, ny=self.ny,
                                    base_h=self.spin_h.value())
            self.view.set_convection(self.H)
            self.log(
                f"Re-applied heatsinks from "
                f"{os.path.basename(self.heatsink_path)} to the "
                f"{self.nx}x{self.ny} grid."
            )
        else:
            self.H = np.full((self.ny, self.nx), self.spin_h.value())
        self.view.set_convection(self.H)

    def show_layout(self):
        """
        Shows the pre-run state: ambient field plus component outlines.

        The view is now a function of state rather than of the action that led
        here, so this is reachable from anywhere instead of being buried inside
        the table-edit handler.
        """
        u_init = np.full((self.ny, self.nx), self.spin_tamb.value())
        # The layout view holds no result, so any probe reading on it would be
        # stale (FIX M-4).
        self.clear_probes(redraw=False)
        self._has_result = False
        self._sim_phase = 'idle'
        self._sim_tmax = None
        self.view.set_temperature(u_init)
        self.view.set_components(self.components_list)
        self.view.set_power(self.Q)
        self.chk_layers['components'].setChecked(True)
        self._sync_layer_controls()
        self.view.set_title(self._compose_title())

    def on_table_edit(self):
        """Called when any cell in the table is edited."""
        self.sync_data_from_table()
        # Cached, not recomputed per redraw, and reusable by the layer views.
        self.Q = calculate_q_matrix(self.components_list, self.nx, self.ny)
        self.show_layout()
        self.log("Component mapping updated.")

    # (column index, component key, value used when the cell is left empty)
    TABLE_NUMERIC_COLUMNS = (
        (1, 'Power_Watts', 0.0),
        (2, 'Center_X_mm', None),
        (3, 'Center_Y_mm', None),
        (4, 'Width_mm', None),
        (5, 'Length_mm', None),
    )
    BAD_CELL_BRUSH = QBrush(QColor(150, 40, 40))
    BAD_CELL_TIP = "Not a number — this component is excluded from the simulation."

    def _mark_cell(self, row, col, valid):
        """Flags a cell that failed to parse. FIX H-6."""
        item = self.table.item(row, col)
        if item is None:
            return
        item.setBackground(QBrush() if valid else self.BAD_CELL_BRUSH)
        item.setToolTip("" if valid else self.BAD_CELL_TIP)

    def sync_data_from_table(self):
        """
        Rebuilds components_list from the table.

        FIX H-6: a typo used to drop the component from the simulation with no
        signal whatsoever — `except (ValueError, AttributeError): continue`.
        The row stayed visible while its power quietly left the model. Bad
        cells are now painted red, given a tooltip, and named in the console;
        the component is still excluded, because its real value is unknown, but
        the exclusion is no longer silent.
        """
        # setBackground writes item data and would re-enter cellChanged.
        self.table.blockSignals(True)
        try:
            self.components_list = []
            problems = []
            for row in range(self.table.rowCount()):
                item = self.table.item(row, 0)
                designator = item.text().strip() if item is not None else ""
                if not designator:
                    designator = f"ROW{row + 1}"

                comp = {'Designator': designator}
                row_valid = True
                for col, key, empty_default in self.TABLE_NUMERIC_COLUMNS:
                    cell = self.table.item(row, col)
                    text = cell.text().strip() if cell is not None else ""
                    value = None
                    if not text and empty_default is not None:
                        value = empty_default
                    elif text:
                        try:
                            value = float(text.replace(",", "."))
                        except ValueError:
                            value = None
                    if value is None:
                        row_valid = False
                        header = self.table.horizontalHeaderItem(col)
                        problems.append(
                            f"{designator}/{header.text() if header else col}"
                            f"={text or '(empty)'!s}"
                        )
                    self._mark_cell(row, col, value is not None)
                    comp[key] = value

                if row_valid:
                    self.components_list.append(comp)

            if problems:
                self.log(
                    f"WARNING: {len(problems)} invalid cell(s) — excluded from "
                    f"the simulation: {', '.join(problems[:6])}"
                    + (" ..." if len(problems) > 6 else "")
                )
        finally:
            self.table.blockSignals(False)

    def add_component_row(self):
        row = self.table.rowCount()
        # FIX M-3: each setItem emits cellChanged, and every one of those used
        # to trigger a full figure rebuild with a fresh colorbar - six for the
        # cells plus the explicit call below, seven redraws for one click.
        self.table.blockSignals(True)
        try:
            self.table.insertRow(row)
            defaults = (f"NEW{row}", "0.0", "50.0", "50.0", "10.0", "10.0")
            for col, value in enumerate(defaults):
                self.table.setItem(row, col, QTableWidgetItem(value))
        finally:
            self.table.blockSignals(False)
        self.on_table_edit()

    def delete_component_row(self):
        current_row = self.table.currentRow()
        if current_row >= 0:
            self.table.removeRow(current_row)
            self.on_table_edit()

    def save_data_to_csv(self):
        import csv
        path, _ = QFileDialog.getSaveFileName(self, "Save Components CSV", "", "CSV Files (*.csv)")
        if path:
            self.sync_data_from_table()
            # FIX H-2: utf-8-sig round-trips through both our reader and Excel.
            with open(path, mode='w', newline='', encoding='utf-8-sig') as file:
                writer = csv.DictWriter(file, fieldnames=["Designator", "Center_X_mm", "Center_Y_mm", "Width_mm", "Length_mm", "Power_Watts"])
                writer.writeheader()
                for comp in self.components_list:
                    writer.writerow(comp)
            self.log(f"Data saved to {path}")

    def start_simulation(self, mode):
        if not self.components_list:
            self.log("Error: No components loaded or defined.")
            return

        # Prepare Q_matrix and K_matrix from current state
        self.Q = calculate_q_matrix(self.components_list, self.nx, self.ny)
        Q = self.Q

        if self.K_matrix is not None:
            K = self.K_matrix
            self.log("Using loaded Gerber for thermal conductivity mapping.")
        else:
            self.log("No Gerber loaded, calculating effective uniform conductivity.")
            # Fallback to uniform K_matrix using k_eff
            k_eff = config.calculate_k_eff(config.BOARD_LAYERS, config.COPPER_OZ, config.K_FR4)
            K = np.full((self.ny, self.nx), k_eff)

        # FIX L-7: read the environment once, here on the GUI thread, and hand
        # it to the job. config.T_amb is no longer written from here, so the
        # worker cannot observe it changing mid-solve.
        t_amb = self.spin_tamb.value()
        h_base = self.spin_h.value()
        emissivity = self.spin_emissivity.value()

        # Reset the field to ambient; the base layer switches back to
        # temperature if the user was looking at conductivity.
        u_init = np.full((self.ny, self.nx), t_amb)
        self.clear_probes(redraw=False)
        self._has_result = True
        self._sim_phase = 'running'
        self._sim_mode = mode
        self._sim_t = 0.0
        self._sim_tmax = None
        # A result is about to arrive: put the base back on temperature if the
        # user was inspecting copper or power.
        self.radio_base['temperature'].setChecked(True)
        self.view.set_temperature(u_init)
        self.view.set_base_kind('temperature')
        self.view.reset_clim()
        self.view.set_clim_lock(self.chk_lock_scale.isChecked(), floor=t_amb)
        self._sync_layer_controls()
        self.view.set_title(self._compose_title())

        # Setup Threading
        self.solver_thread = QThread()
        self.worker = SolverWorker(self.nx, self.ny, Q,
                                   config.SELECTED_MATERIAL.split(" ")[0],
                                   self.H, config.BOARD_LAYERS, config.COPPER_OZ,
                                   mode=mode, t_final=self.spin_tfinal.value(),
                                   K_matrix=K,
                                   T_amb=t_amb, h_ambient=h_base,
                                   emissivity=emissivity)
        
        self.worker.moveToThread(self.solver_thread)
        
        # Connect signals
        self.solver_thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.on_simulation_finished)
        self.worker.finished.connect(self.solver_thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.solver_thread.finished.connect(self.solver_thread.deleteLater)
        
        self.worker.progress.connect(self.update_plots)
        self.worker.log.connect(self.log)

        # FIX H-5: lock every control, not just the two run buttons.
        self._set_busy(True, can_stop=True)

        self.solver_thread.start()

    def stop_simulation(self):
        if self.worker is None:
            return
        try:
            self.worker.stop()
        except RuntimeError:
            # The worker's C++ side was already deleted; nothing left to stop.
            self.worker = None
            return
        self.btn_stop.setEnabled(False)
        self.log("Stopping simulation...")

    # ----- Shutdown -----

    SHUTDOWN_TIMEOUT_MS = 5000

    def closeEvent(self, event):
        """
        Shuts the solver thread down before the window goes away.

        FIX C-5: without this, closing the window during a run left the QThread
        alive. Qt then destroyed the QThread wrapper while the thread was still
        executing, which prints "QThread: Destroyed while thread is still
        running" and aborts the process. A steady-state run made it worse: it
        polled no cancellation flag at all, so there was nothing to stop.

        worker.stop() is what actually breaks the compute loop - quit() only
        unwinds the thread's event loop once run() has returned. quit() is
        posted to the worker thread's own queue, so it is picked up as soon as
        run() finishes, without needing the (blocked) GUI event loop.
        """
        try:
            self._shutdown_solver_thread()
            # The Gerber worker has no cancellation point, so this can only
            # wait it out; PyGerber rasterisation is seconds, not minutes.
            self._shutdown_thread(self.gerber_thread, "Gerber loader")
        finally:
            # Always detach, on every exit path: the handler is attached to
            # module-level loggers that outlive this window, and a stale handler
            # would keep writing into a destroyed QTextEdit.
            self._detach_logging()
        event.accept()

    def _shutdown_solver_thread(self):
        try:
            if self.worker is not None and self.solver_thread is not None \
                    and self.solver_thread.isRunning():
                self.log("Closing: stopping the running simulation...")
                self.worker.stop()
        except RuntimeError:
            pass
        self._shutdown_thread(self.solver_thread, "Solver")

    def _shutdown_thread(self, thread, label):
        """Waits for one worker thread, terminating it only as a last resort."""
        if thread is None:
            return
        try:
            if not thread.isRunning():
                return
        except RuntimeError:
            return  # wrapper already deleted, nothing to wait for

        thread.quit()
        if not thread.wait(self.SHUTDOWN_TIMEOUT_MS):
            self.logger.warning(
                "%s thread did not stop within %d ms; terminating it.",
                label, self.SHUTDOWN_TIMEOUT_MS,
            )
            thread.terminate()
            thread.wait(2000)

    def _detach_logging(self):
        """Stops feeding records to a console widget that is going away."""
        for name in self.LOGGED_MODULES:
            logging.getLogger(name).removeHandler(self.handler)

    def update_plots(self, u, t, max_t):
        # FIX M-4: probes read whatever is on screen, so they track a transient
        # instead of reporting the previous run's numbers.
        self._sim_t = t
        self._sim_tmax = max_t
        self.view.set_temperature(u)
        self.view.set_title(self._compose_title())
        self._refresh_probe_labels()
        self.view.request_draw()

    def on_simulation_finished(self, u_final):
        self._set_busy(False)

        # The references are deliberately NOT cleared here. This slot runs
        # inside the emission of worker.finished, and self.worker holds the last
        # Python reference to the worker: dropping it destroys the C++ object
        # while Qt is still walking the remaining connections (quit,
        # deleteLater), which crashes the process with STATUS_STACK_BUFFER_OVERRUN.
        # Stale references are handled by the RuntimeError guards in
        # stop_simulation() and _shutdown_solver_thread() instead.
        if u_final is not None:
            self.u_final = u_final
            self._sim_phase = 'done'
            # Probes are NOT cleared here any more: they are cleared at the
            # start of a run, and any placed during a transient should settle
            # onto the final values rather than vanish at the finish (FIX M-4).
            self.update_plots(u_final, self.spin_tfinal.value(), np.max(u_final))
            self.log("Simulation Result Rendered.")
        else:
            self._sim_phase = 'idle'
            self._has_result = False
        self._sync_layer_controls()

    # ----- Interactive Virtual Probes -----

    def on_canvas_click(self, event):
        """
        Handles mouse clicks on the Matplotlib canvas.
        Left-click: place a virtual probe (max 10).
        Right-click: remove the nearest probe.
        """
        # Probing during a transient is allowed and useful (FIX M-4): both the
        # click and the worker's progress signal are handled on the GUI thread,
        # so they serialise. What must not happen concurrently is a figure
        # rebuild, and the table that triggers one is locked while busy (H-5).
        if (event.inaxes != self.ax
                or not self._has_result
                or self.view.field_kind() != 'temperature'
                or self.view.field() is None):
            return

        x_mm, y_mm = event.xdata, event.ydata
        if x_mm is None or y_mm is None:
            return

        if event.button == 1:  # Left click -> place probe
            if len(self.probes) >= self.MAX_PROBES:
                self.log(f"Max probes ({self.MAX_PROBES}) reached. Right-click to remove or Clear Probes.")
                return
            self._place_probe(x_mm, y_mm)

        elif event.button == 3:  # Right click -> remove nearest probe
            self._remove_nearest_probe(x_mm, y_mm)

    def _probe_indices(self, x_mm, y_mm):
        """
        Grid cell containing a physical point.

        FIX M-4: this used round(). With origin='lower' and extent=[0, W, 0, H]
        cell i covers [i*dx, (i+1)*dx], so its centre sits at (i+0.5)*dx.
        Rounding therefore picked the neighbouring cell for every point in the
        lower half of a cell - a systematic half-cell (0.25 mm) offset. floor,
        which is what int() does for the non-negative side, picks the cell the
        point is actually inside.
        """
        dx_mm = config.dx * 1000.0
        # Clicks arrive in CAD millimetres (the axes carry the CAD frame), so
        # the board origin comes off before the index is taken.
        ix = int((x_mm - config.BOARD_ORIGIN_X_MM) / dx_mm)
        iy = int((y_mm - config.BOARD_ORIGIN_Y_MM) / dx_mm)
        return (max(0, min(self.ny - 1, iy)), max(0, min(self.nx - 1, ix)))

    def _refresh_probe_labels(self):
        """Re-reads every probe from the field currently displayed (FIX M-4)."""
        # Only the temperature field is in degrees; refusing to relabel from a
        # conductivity or power base is what keeps the unit on the label true.
        if self.view.field_kind() != 'temperature':
            return
        field = self.view.field()
        if field is None:
            return
        for (px, py), (_marker, label) in zip(self.probes, self.probe_artists):
            iy, ix = self._probe_indices(px, py)
            try:
                label.set_text(f"{field[iy, ix]:.1f} °C")
            except (ValueError, RuntimeError, IndexError):
                pass  # artist detached by a figure rebuild

    def _place_probe(self, x_mm, y_mm):
        """Places a probe marker + temperature annotation at (x_mm, y_mm)."""
        iy, ix = self._probe_indices(x_mm, y_mm)
        temperature = self.view.field()[iy, ix]

        # Draw marker and label
        marker, = self.ax.plot(
            x_mm, y_mm, 'x', color='cyan', markersize=10,
            markeredgewidth=2, zorder=BoardView.Z_PROBES
        )
        label = self.ax.annotate(
            f"{temperature:.1f} \u00b0C",
            xy=(x_mm, y_mm),
            xytext=(8, 8), textcoords='offset points',
            fontsize=9, fontweight='bold',
            color='white',
            bbox=dict(boxstyle='round,pad=0.3', fc='black', alpha=0.7),
            zorder=BoardView.Z_PROBES + 1
        )

        self.probes.append((x_mm, y_mm))
        self.probe_artists.append((marker, label))
        self._publish_probe_artists()

        probe_num = len(self.probes)
        self.log(f"Probe #{probe_num}: ({x_mm:.1f}, {y_mm:.1f}) mm = {temperature:.1f} \u00b0C")

    def _remove_nearest_probe(self, x_mm, y_mm):
        """Removes the probe closest to the click position."""
        if not self.probes:
            return

        # Find nearest probe
        dists = [(px - x_mm)**2 + (py - y_mm)**2 for px, py in self.probes]
        idx = int(np.argmin(dists))

        # Remove artists from plot
        marker, label = self.probe_artists[idx]
        try:
            marker.remove()
            label.remove()
        except (ValueError, NotImplementedError):
            pass

        removed = self.probes.pop(idx)
        self.probe_artists.pop(idx)
        self._publish_probe_artists()
        self.log(f"Removed probe at ({removed[0]:.1f}, {removed[1]:.1f}) mm.")

    def clear_probes(self, redraw=True):
        """Removes all probe markers from the canvas."""
        for marker, label in self.probe_artists:
            try:
                marker.remove()
                label.remove()
            except (ValueError, NotImplementedError):
                pass  # Artist already removed
        self.probes.clear()
        self.probe_artists.clear()
        self._publish_probe_artists()
        if redraw:
            self.view.request_draw()

    # ----- Save Result Image -----

    def save_result_image(self):
        """Saves the current figure (with probes) to a PNG/JPEG file."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Result Image", "",
            "PNG Image (*.png);;JPEG Image (*.jpg *.jpeg);;All Files (*)"
        )
        if not path:
            return
        # Redraws are deferred, so make sure the figure matches the state
        # before it is rendered to disk.
        self.view.flush_pending()
        self.figure.savefig(path, dpi=200, bbox_inches='tight')
        self.log(f"Result image saved to: {path}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
