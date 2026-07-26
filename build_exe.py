import PyInstaller.__main__
import os

import config

# Формируем аргументы для PyInstaller (используем ';' как разделитель для Windows)
args = [
    'main.py',
    f'--name=HFDM_Solver_v{config.VERSION}',
    '--onefile',
    '--windowed',
    '--clean',
    '--hidden-import=pygerber'
]

# Запускаем программную сборку
PyInstaller.__main__.run(args)
