# run as: python py2exe.py build
from cx_Freeze import setup, Executable
import os
import shutil

setup(
    name='PUL BLAST',
    version='1',
    options={'build_exe': {"includes": ["tkinter", "Bio", "xml"]}},
    executables=[Executable('gui.py',
    targetName='pulBlast',
    copyright='David Miškić')]
)

if os.path.isdir("build/exe.linux-x86_64-3.6/lib/Tkinter"):
    try:
        os.rename("build/exe.linux-x86_64-3.6/lib/Tkinter", "build/exe.linux-x86_64-3.6/lib/tkinter")
    except Exception as e:
        print(e)
        print("Renaming Tkinter folder failed. Check manually.")

try:
    shutil.copytree(src="javascript", dst="build/exe.linux-x86_64-3.6/javascript")
    shutil.copyfile(src="clankiDB.phr", dst="build/exe.linux-x86_64-3.6/clankiDB.phr")
    shutil.copyfile(src="clankiDB.pin", dst="build/exe.linux-x86_64-3.6/clankiDB.pin")
    shutil.copyfile(src="clankiDB.psq", dst="build/exe.linux-x86_64-3.6/clankiDB.psq")
    shutil.copyfile(src="modularityDescription.txt", dst="build/exe.linux-x86_64-3.6/modularityDescription.txt")
except Exception as e:
    print(e)
    print("Copying additional files failed. Check manually.")