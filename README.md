[![Twitter Follow](https://img.shields.io/twitter/follow/plumed_org.svg?style=social&label=Follow)](https://twitter.com/mdpoleto)

# RAIJIN : Electric field analysis for molecular simulations

Raijin is a python algorithm that employs MDanalysis engine to calculate Electric Field at any point inside
the simulation box throughout MD trajectories. Raijin also includes a PyMOL plugin to visualize electric
field vectors in your simulation box.

Required packages:

* MDanalysis >= 1.0.0
* Python     >= 3.x
* Numpy      >= 1.2.x


Usage and Examples
------------------------------

We have included examples of all 3 MODES implemented in Raijin.

* ATOM mode:

* BOND mode:

* COORDINATE mode:


Installation instructions
------------------------------

First, make sure you have all required packages installed. For MDanalysis installation procedures, [click here](https://www.mdanalysis.org/pages/installation_quick_start/).

After, just clone this repository into a folder of your choice:

    git clone https://github.com/mdpoleto/raijin.git

To use Raijin easily, copy the directory pathway to Raijin folder and include an alias in your ~/.bashrc:

    alias raijin="python /path/to/the/cloned/repository/raijin.py"
