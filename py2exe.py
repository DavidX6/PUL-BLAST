# run as: python py2exe.py build
from cx_Freeze import setup, Executable
import os
import shutil

setup(
    name='PUL BLAST',
    version='1',
    options={'build_exe': {"includes": ["tkinter", "Bio", "xml"]}},
    executables=[Executable('gui.py',
                            targetName='pulBlast.exe',
                            copyright='David Miškić',
                            base='Win32GUI')]
)
if os.path.isdir("build\exe.win32-3.8\lib\Tkinter"):
    try:
        os.rename("build/exe.win32-3.8/lib/Tkinter", "build/exe.win32-3.8/lib/tkinter")
    except Exception as e:
        print(e)
        print("Renaming Tkinter folder failed. Check manually.")

try: shutil.copytree(src="javascript", dst="build\exe.win32-3.8\javascript")
except Exception as e: print(e)
try: shutil.copyfile(src="clankiDB.phr", dst="build\exe.win32-3.8\clankiDB.phr")
except Exception as e: print(e)
try: shutil.copyfile(src="clankiDB.pin", dst="build\exe.win32-3.8\clankiDB.pin")
except Exception as e: print(e)
try: shutil.copyfile(src="clankiDB.psq", dst="build\exe.win32-3.8\clankiDB.psq")
except Exception as e: print(e)
try: shutil.copyfile(src="modularityDescription.txt", dst="build\exe.win32-3.8\modularityDescription.txt")
except Exception as e: print(e)