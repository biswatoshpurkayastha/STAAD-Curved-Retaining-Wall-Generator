Curved Retaining Wall Generator for STAAD.Pro
Overview
Curved retaining walls are difficult to model by hand in STAAD.Pro. The engineer has to:

Build a curved wall surface from plate (shell) elements
Apply water pressure and earth pressure, both of which vary with depth
Represent the soil behind or below the wall, which can push on the wall but cannot pull it
Doing this manually is slow and easy to get wrong, especially when the radius, height or soil data changes.

This tool automates the whole process. You enter the wall geometry, soil properties and water level in a desktop interface (GUI), preview the wall in 3D, and then generate the model directly in STAAD.Pro through OpenStaadPy.

<img src="screenshot_gui.png" alt="Curved retaining wall GUI" width="800">

Key Features
Desktop interface built with PySide6
Interactive 3D preview of the curved wall mesh using PyVista, before anything is sent to STAAD.Pro
Curved wall profile generated as STAAD.Pro plate (shell) elements
Water and earth pressure applied as loads that vary with depth
Soil–structure interaction modelled with compression-only members, so the soil supports the wall in compression only
Direct model creation in STAAD.Pro through the OpenStaadPy library
Workflow
Open STAAD.Pro and create or open a blank model.
Run the Python script.
Enter the wall geometry, soil properties and water level in the GUI.
Check the wall in the 3D PyVista preview.
Click Generate to create nodes, plates, supports and loads in STAAD.Pro.
Run the analysis in STAAD.Pro and review the results.
Requirements
Windows
STAAD.Pro with OpenSTAAD available (Python bitness must match STAAD.Pro, normally 64-bit)
Python 3.x
Python libraries listed in requirements.txt
Installation
Clone or download this repository, then open Command Prompt or PowerShell in the project folder and run:

python -m pip install -r requirements.txt
If you have more than one Python installation, use the full path to the Python you want to run, for example:

& "C:\Program Files\Python313\python.exe" -m pip install --user -r requirements.txt
How to Run
Start STAAD.Pro and open a model.
Run:
python retaining_wall_gui.py
Troubleshooting
Error	Fix
ModuleNotFoundError: No module named 'PySide6'	python -m pip install PySide6
ModuleNotFoundError: No module named 'pyvista'	python -m pip install pyvista pyvistaqt
Script cannot connect to STAAD.Pro	Make sure STAAD.Pro is open with a model loaded before you run the script
Package installed but error persists	Install with the same Python that runs the script: <path-to-python.exe> -m pip install ...
Libraries Used
Purpose	Library
STAAD.Pro integration	OpenStaadPy
User interface	PySide6
3D visualisation	PyVista, PyVistaQt
Numerical computation	NumPy
Disclaimer
This tool is intended to help with modelling. Engineers must check the generated geometry, loads, supports and soil parameters before using the results for design.

Author
Biswatosh Purkayastha
