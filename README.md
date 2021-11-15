[![Twitter Follow](https://img.shields.io/twitter/follow/plumed_org.svg?style=social&label=Follow)](https://twitter.com/mdpoleto)

# RAIJIN : Electric field analysis for molecular simulations

## What is ***Raijin***?
Raijin is a python algorithm that employs MDanalysis engine to calculate Electric Field at any point inside
the simulation box throughout MD trajectories. Raijin also includes a PyMOL plugin to visualize electric
field vectors in your simulation box.

Required packages:

* MDanalysis >= 1.0.0
* Python     >= 3.x
* Numpy      >= 1.2.x


##Installation instructions
------------------------------

First, make sure you have all required packages installed. For MDanalysis installation procedures, [click here](https://www.mdanalysis.org/pages/installation_quick_start/).

After, just clone this repository into a folder of your choice:

    git clone https://github.com/mdpoleto/raijin.git

To use Raijin easily, copy the directory pathway to Raijin folder and include an alias in your ~/.bashrc:

    alias raijin="python /path/to/the/cloned/repository/raijin.py"

To install the PyMOL plugin, open PyMOL > Plugin Manager and click on "Install New Plugin" tab.
Load the ***Raijin*** plugin and use it via command-line within PyMOL. To usage instructions, read our FAQ.


## Usage and Examples
------------------------------
***Raijin*** calculations are based on parameters that are provided via a configuration file,
which can be obtained via the command:

   raijin -template

***Raijin*** has 3 calculations MODES:

* In ***ATOM*** mode, the coordinate of one atom will be tracked throughout the trajectory to serve as target point.

* In ***BOND*** mode, the midpoint between 2 atoms will be tracked throughout the trajectory to serve as target
point. In ***BOND*** mode, the bond axis is used to calculate electric field alignment. By default, the bond axis is
define as **selbond1 ---> selbond2**.

* In ***COORDINATE*** mode, a list of [X,Y,Z] coordinates will serve as target point in all trajectory frames.
