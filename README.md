The Arbitrary Composition Condensation Sequence Calculator (ArCCoS) is an equilibrium condensation sequence calculator. It calculates the distribution of 23 elements between ~400 gas phases and ~100 solid phases of interest to planetary and exoplanetary science. See https://arxiv.org/abs/1604.08309 more information

ArCCoS is released under the GNU GPL v2 or newer

Source code: https://github.com/CaymanUnterborn/ArCCoS

Authors (as of 2016, listed alphabetically by first name): Dr. Cayman Unterborn (unterborn.1@osu.edu)

Basic Structure: The code is split into 5 files: Input (input.py) Main code (arccos.py) Mass balance function (condensation/fun.py) Database calculations & entropy/enthalpy calculator(condensation/get_data.py) Output (condensation/write.py)

as well as a database of thermodynamic values for the gasses/solids (Data/).

Following the prompts in input.py will show the current input variables including: Solar composition model and total pressure. Future versions will include database manipulation as well as more complex thermodynamics (e.g. solid solutions).

After changing the input file, to run ArCCoS, simply type: python input.py and the code will run until 100% of all refractory elements are condensed. At completion, matplotlib will produce a plot showing the refractory element condensation (% element condensed vs. temperature). Further output files including a text version of the appearance and disappearance temperature of each solid (output/sequence/) and individual element distributions as a function of temperature (output/abundance/).

Requirements

Python 2.7.x or Python 3.4+ Python modules: NumPy, SciPy, Matplotlib

Install under Ubuntu

Install using apt by opening a terminal window and entering sudo apt-get install python python-scipy python-numpy python-matplotlib
Go to the Burnman examples directory and type: python example_beginner.py Figures should show up, indicating that it is working.
Install on a Mac

get Xcode
If you don't have Python yet, download it (for free) from python.org/download . Make sure to use either Python 2.7 or Python 3.4+. To check your version of python, type the following in a terminal: python --version
Install the latest Numpy version: http://sourceforge.net/projects/numpy/files/NumPy/
Install the latest Scipy at http://sourceforge.net/projects/scipy/files/
Install the latest Matplotlib from http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.1.1/
Go to the Burnman examples directory and type: python example_beginner.py Figures should show up, indicating that it is working.
Python can also be downloaded in an ready to use fashion by using the Canopy distribution of python located at: https://www.enthought.com/products/canopy/ with numpy, scipy and matplotlib installs being accessed through the Package Manager.

Problems you might run into:

Installing numpy/scipy/matplotlib for a different python version than the one on your computer

Having matplotlib for 32-bit instead of 64-bit (for me this got fixed by installing the very latest version). This will give you the error no matching architecture in universal wrapper. You can check if your python distribution is 32 or 64 bit with the following lines:

python 
>>> import platform
>>> print platform.architecture()
Install under Windows

To get Python 2.7.x (for example) running under Windows:

Download Python from http://www.python.org/ and install the version at C:\Python27\; the 32-bit version is recommended
Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy, download "numpy-MKL-1.6.2.win32-py2.7.exe" and install
Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy, download "scipy-0.10.1.win32-py2.7.exe" and install
Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#matplotlib, download "matplotlib-1.1.1.win32-py2.7.exe" and install
Open Python Shell (IDLE Python GUI)
File -- Open -- find one of the example files
Run the module (or press F5)
