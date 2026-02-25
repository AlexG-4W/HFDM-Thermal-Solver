import sys
import logging
import numpy as np
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QPushButton, QFileDialog, QGroupBox, 
                             QLabel, QDoubleSpinBox, QTextEdit, QFrame, QSplitter,
                             QTableWidget, QTableWidgetItem, QTabWidget)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QObject
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

from data_loader import (load_components_dict, load_heatsinks, 
                         load_components_list, calculate_q_matrix)
from solver import PCBSolver
from visualization import plot_heatmap, update_heatmap, plot_probe_history
import config

# --- Logging Handler for GUI ---
class QTextEditHandler(logging.Handler):
    def __init__(self, text_edit):
        super().__init__()
        self.text_edit = text_edit

    def emit(self, record):
        msg = self.format(record)
        self.text_edit.append(msg)

# --- Worker for Background Calculations ---
class SolverWorker(QObject):
    finished = pyqtSignal(object)  # Emits final u matrix
    progress = pyqtSignal(object, float, float) # Emits (u, current_time, max_temp)
    log = pyqtSignal(str)

    def __init__(self, nx, ny, Q, material, h_matrix, layers, copper_oz, mode='steady', t_final=100.0, K_matrix=None):
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
                               K_matrix=self.K_matrix)
            
            self.log.emit(f"Starting {self.mode} simulation...")
            self.log.emit(f"Average k: {np.mean(solver.K):.4f} W/mK")

            if self.mode == 'steady':
                u_final, iterations = solver.solve_steady_state()
                self.log.emit(f"Steady-state converged in {iterations} iterations.")
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

# --- Main Window ---
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("HFDM Thermal Solver - Desktop")
        self.resize(1200, 800)

        self.comp_dict = None
        self.components_list = []
        # Calculate nx, ny based on board size (100x100mm) and dx (0.5mm)
        self.nx = int(100.0 / (config.dx * 1000))
        self.ny = int(100.0 / (config.dx * 1000))
        self.H = None
        self.K_matrix = None
        self.solver_thread = None
        self.worker = None

        self.init_ui()
        self.setup_logging()

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
        
        # Gerber Input
        self.btn_load_gerber = QPushButton("Load Top Copper (Gerber)")
        self.btn_load_gerber.clicked.connect(self.load_gerber)
        data_layout.addWidget(self.btn_load_gerber)
        
        self.lbl_gerber_status = QLabel("Gerber: None loaded")
        data_layout.addWidget(self.lbl_gerber_status)
        
        self.btn_view_topo = QPushButton("View Topology")
        self.btn_view_topo.clicked.connect(self.view_topology)
        self.btn_view_topo.setEnabled(False)
        data_layout.addWidget(self.btn_view_topo)

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
        self.table.setHorizontalHeaderLabels(["Designator", "Power [W]", "X [mm]", "Y [mm]", "W [mm]", "H [mm]"])
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
        
        self.figure = Figure(figsize=(8, 6))
        self.canvas = FigureCanvas(self.figure)
        right_layout.addWidget(self.canvas)
        
        self.ax = self.figure.add_subplot(111)
        self.im = None
        
        splitter.addWidget(right_panel)
        
        main_layout.addWidget(splitter)

    def setup_logging(self):
        self.handler = QTextEditHandler(self.log_console)
        self.handler.setFormatter(logging.Formatter('%(asctime)s - %(message)s', '%H:%M:%S'))
        self.logger = logging.getLogger("HFDM_GUI")
        self.logger.addHandler(self.handler)
        self.logger.setLevel(logging.INFO)
        self.logger.info("HFDM GUI Initialized.")

    def log(self, message):
        self.logger.info(message)

    def load_data(self):
        # In a real app, we'd use a file dialog. For now, we use the default loader logic.
        try:
            self.components_list = load_components_list()
            # Ensure nx, ny are fresh based on current config
            self.nx = int(100.0 / (config.dx * 1000))
            self.ny = int(100.0 / (config.dx * 1000))
            self.H = load_heatsinks(nx=self.nx, ny=self.ny)
            
            self.populate_table()
            
            self.lbl_status.setText(f"Loaded: {self.nx}x{self.ny} grid")
            self.log(f"Successfully loaded {len(self.components_list)} components and heatsinks.")
            self.btn_save_csv.setEnabled(True)
            
            # Show initial empty board or placeholder
            self.on_table_edit() # This triggers initial rendering
            
        except Exception as e:
            self.log(f"Load Error: {str(e)}")

    def load_gerber(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load Top Copper Gerber", "", "Gerber Files (*.gbr *.gtl);;All Files (*)")
        if path:
            try:
                from data_loader import load_gerber_to_k_matrix
                import os
                self.log("Loading Gerber file. This may take a moment...")
                QApplication.processEvents()
                
                # Assume 100x100 board
                self.K_matrix = load_gerber_to_k_matrix(path, 100.0, 100.0, config.dx, config.K_FR4, config.K_CU)
                filename = os.path.basename(path)
                self.lbl_gerber_status.setText(f"Gerber: {filename}")
                self.btn_view_topo.setEnabled(True)
                self.log(f"Successfully loaded Gerber: {filename}")
            except Exception as e:
                self.log(f"Gerber Load Error: {str(e)}")

    def view_topology(self):
        if self.K_matrix is not None:
            self.figure.clear()
            self.ax = self.figure.add_subplot(111)
            # Visualize K_matrix. High K = copper, low K = FR-4
            im = self.ax.imshow(self.K_matrix, cmap='copper', origin='lower', extent=[0, 100.0, 0, 100.0])
            self.figure.colorbar(im, ax=self.ax, label='Thermal Conductivity [W/mK]')
            self.ax.set_title("Copper Topology (Thermal Conductivity)")
            self.ax.set_xlabel('X [mm]')
            self.ax.set_ylabel('Y [mm]')
            self.im = im
            self.canvas.draw()
            self.log("Displaying copper topology.")

    def populate_table(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.components_list))
        for i, comp in enumerate(self.components_list):
            self.table.setItem(i, 0, QTableWidgetItem(str(comp['Designator'])))
            self.table.setItem(i, 1, QTableWidgetItem(str(comp['Power_Watts'])))
            self.table.setItem(i, 2, QTableWidgetItem(str(comp['Center_X_mm'])))
            self.table.setItem(i, 3, QTableWidgetItem(str(comp['Center_Y_mm'])))
            self.table.setItem(i, 4, QTableWidgetItem(str(comp['Width_mm'])))
            self.table.setItem(i, 5, QTableWidgetItem(str(comp['Height_mm'])))
        self.table.blockSignals(False)

    def on_table_edit(self):
        """Called when any cell in the table is edited."""
        self.sync_data_from_table()
        
        # Calculate Q_matrix
        Q = calculate_q_matrix(self.components_list, self.nx, self.ny)
        
        # Update visualization with initial temperature
        u_init = np.full((self.ny, self.nx), self.spin_tamb.value())
        
        # IMPORTANT: Use clear_fig=True to prevent colorbar accumulation
        self.im = plot_heatmap(u_init, 100.0, 100.0, self.figure, self.ax, clear_fig=True)
        # Update self.ax because fig.clear() might have invalidated it
        self.ax = self.figure.axes[0]
        self.canvas.draw()
        self.log("Component mapping updated.")

    def sync_data_from_table(self):
        self.components_list = []
        for i in range(self.table.rowCount()):
            try:
                comp = {
                    'Designator': self.table.item(i, 0).text() if self.table.item(i, 0) else f"C{i}",
                    'Power_Watts': float(self.table.item(i, 1).text()) if self.table.item(i, 1) else 0.0,
                    'Center_X_mm': float(self.table.item(i, 2).text()) if self.table.item(i, 2) else 0.0,
                    'Center_Y_mm': float(self.table.item(i, 3).text()) if self.table.item(i, 3) else 0.0,
                    'Width_mm': float(self.table.item(i, 4).text()) if self.table.item(i, 4) else 1.0,
                    'Height_mm': float(self.table.item(i, 5).text()) if self.table.item(i, 5) else 1.0
                }
                self.components_list.append(comp)
            except (ValueError, AttributeError):
                continue

    def add_component_row(self):
        row = self.table.rowCount()
        self.table.insertRow(row)
        self.table.setItem(row, 0, QTableWidgetItem(f"NEW{row}"))
        for j in range(1, 6):
            self.table.setItem(row, j, QTableWidgetItem("0.0"))
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
            with open(path, mode='w', newline='') as file:
                writer = csv.DictWriter(file, fieldnames=["Designator", "Center_X_mm", "Center_Y_mm", "Width_mm", "Height_mm", "Power_Watts"])
                writer.writeheader()
                for comp in self.components_list:
                    writer.writerow(comp)
            self.log(f"Data saved to {path}")

    def start_simulation(self, mode):
        if not self.components_list:
            self.log("Error: No components loaded or defined.")
            return

        # Prepare Q_matrix and K_matrix from current state
        Q = calculate_q_matrix(self.components_list, self.nx, self.ny)
        
        if self.K_matrix is not None:
            K = self.K_matrix
            self.log("Using loaded Gerber for thermal conductivity mapping.")
        else:
            self.log("No Gerber loaded, calculating effective uniform conductivity.")
            # Fallback to uniform K_matrix using k_eff
            k_eff = config.calculate_k_eff(config.BOARD_LAYERS, config.COPPER_OZ, config.K_FR4)
            K = np.full((self.ny, self.nx), k_eff)

        # Update config with GUI values
        config.T_amb = self.spin_tamb.value()

        # Reset plot to initial temperature to clear topography colorbar
        u_init = np.full((self.ny, self.nx), self.spin_tamb.value())
        self.im = plot_heatmap(u_init, 100.0, 100.0, self.figure, self.ax, clear_fig=True)
        self.ax = self.figure.axes[0]
        self.canvas.draw()

        # Setup Threading
        self.solver_thread = QThread()
        self.worker = SolverWorker(self.nx, self.ny, Q, 
                                   config.SELECTED_MATERIAL.split(" ")[0], 
                                   self.H, config.BOARD_LAYERS, config.COPPER_OZ,
                                   mode=mode, t_final=self.spin_tfinal.value(),
                                   K_matrix=K)
        
        self.worker.moveToThread(self.solver_thread)
        
        # Connect signals
        self.solver_thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.on_simulation_finished)
        self.worker.finished.connect(self.solver_thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.solver_thread.finished.connect(self.solver_thread.deleteLater)
        
        self.worker.progress.connect(self.update_plots)
        self.worker.log.connect(self.log)

        # UI State
        self.btn_steady.setEnabled(False)
        self.btn_transient.setEnabled(False)
        self.btn_stop.setEnabled(True)

        self.solver_thread.start()

    def stop_simulation(self):
        if self.worker:
            self.worker.stop()
            self.log("Stopping simulation...")

    def update_plots(self, u, t, max_t):
        if self.im is None:
            self.im = plot_heatmap(u, 100.0, 100.0, self.figure, self.ax)
        else:
            update_heatmap(self.im, u, self.ax, f"Time: {t:.1f}s | Max Temp: {max_t:.1f} °C")
        self.canvas.draw()

    def on_simulation_finished(self, u_final):
        self.btn_steady.setEnabled(True)
        self.btn_transient.setEnabled(True)
        self.btn_stop.setEnabled(False)
        
        if u_final is not None:
            self.update_plots(u_final, self.spin_tfinal.value(), np.max(u_final))
            self.log("Simulation Result Rendered.")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
