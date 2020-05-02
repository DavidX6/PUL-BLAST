from cx_Freeze import setup, Executable
import os.path

PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')

include_files = [os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tk86t.dll'),
		 os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tcl86t.dll')]

setup(
    name='PUL BLAST',
    description='A brief description of your program.',
    version='1',
    options={'build_exe': {'include_files': include_files, "includes":["tkinter", "Bio", "xml"]}},
    executables=[Executable('gui.py',
    targetName='pulBlast.exe',
    copyright='Copyright (C) <name> 2018',
    base='Win32GUI')]
)