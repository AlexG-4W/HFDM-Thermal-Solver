import sys
from PyQt6.QtWidgets import QApplication
from gui import MainWindow

def main():
    """
    Entry point for the HFDM Thermal Solver GUI application.
    """
    app = QApplication(sys.argv)
    
    # Optional: Set application style or icon here
    # app.setStyle("Fusion")
    
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
