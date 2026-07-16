import PyInstaller.__main__
import os

# Формируем аргументы для PyInstaller (используем ';' как разделитель для Windows)
args = [
    'main.py',
    '--name=HFDM_Solver_v1.1',
    '--onefile',
    '--windowed',
    '--clean',
    '--hidden-import=pygerber'
]

# Запускаем программную сборку
PyInstaller.__main__.run(args)
