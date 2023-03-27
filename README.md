About the necessary Python 3.6+ compiler:
----

On Windows its possible to run python programms on Visual Studio Code (see   https://code.visualstudio.com/docs/setup/windows).

On Linux S.O. its easier to run python, as on the lastest versions of Ubuntu, python is pre-installed and you can check its version with:
---

$ python --version

If you don't have python on your system its quite easy to install with the following:

$ sudo apt update

$ sudo apt install build-essential

$ sudo apt install python3.8

If you are having problems, there are multiple easy to follow guides like these: https://phoenixnap.com/kb/how-to-install-python-3-ubuntu#:~:text=How%20to%20Install%20Python%203%20on%20Ubuntu%2018.04,Code%20%28Latest%20Version%29%20Step%201%3A%20Update%20Local%20Repositories

https://blog.eldernode.com/install-python-3-ubuntu-20/

------------------------------------------------------------------------------------------

#Molecular-Dynamics related files.
---

AvR_DockingPlot.py - This code is a matplolib script for a very specific analysis and is useless for most, but it could help you if you are not that familiar with matplotlib.


The following was tested on both Windows and Linux S.O.:
---

Analysis_Plots_4.0.py - RMSD, RMSF and Radgyr data analysis. This code needs numpy and matplotlib modules with an python 3.6+ compiler (if you dont usually use Terminal, install Visual Studio Code, which runs on Windows and Linux). 

--------------------------------------
About Analysis_Plots.py execution, You will need some libraries:

https://pypi.org/project/numpy/

https://pypi.org/project/scipy/

https://pypi.org/project/matplotlib/

if it still did not work, try to install

https://pypi.org/project/mpl-toolkits.clifford/


---
The following was tested only on Linux S.O. or on the Windows Subsystem for Linux (WSL):
---

Multidock_Vna.py - This code makes it easy to rerun Vina as many times as you want. Docking uses sampling techniques so trying a few times is better than accepting the first "best model" given by vina.

ASM.py - Molecular Dynamics simulation manager (tested on AMBER18 and AMBERtools18+). Automates the process of performing a molecular dynamics simulation using the AMBER/AMBERtools computational packages.

  Ps: AMBERtools is an open source package of AMBER and both can be found on https://ambermd.org/
  
This programm needs both Python 3.6+ and AMBER/AMBERtools installed on the running machine.
