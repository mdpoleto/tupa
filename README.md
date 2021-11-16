[![Twitter Follow](https://img.shields.io/twitter/follow/plumed_org.svg?style=social&label=Follow)](https://twitter.com/mdpoleto)

# RAIJIN : Electric field analysis for molecular simulations

## What is ***Raijin***?
***Raijin*** is a python algorithm that employs MDanalysis engine to calculate Electric Field at any point inside
the simulation box throughout MD trajectories. ***Raijin*** also includes a PyMOL plugin to visualize electric
field vectors in your simulation box.

Required packages:

* MDanalysis >= 1.0.0
* Python     >= 3.x
* Numpy      >= 1.2.x


## Installation instructions
------------------------------

First, make sure you have all required packages installed. For MDanalysis installation procedures, [click here](https://www.mdanalysis.org/pages/installation_quick_start/).

After, just clone this repository into a folder of your choice:

    git clone https://github.com/mdpoleto/raijin.git

To easily call ***Raijin***, copy the directory pathway to its folder and include an alias in your ~/.bashrc:

    alias raijin="python /path/to/the/cloned/repository/raijin.py"

## Raijin Usage
------------------------------
***Raijin*** calculations are based on parameters that are provided via a configuration file,
which can be obtained via the command:

    raijin -template config.conf


The configuration file usually contains:
    [Elecfield Selection]
    sele_elecfield      = (string)             [default: None]

    [Target Selection]
    mode                = (string)             [default: None]
    selatom             = (string)             [default: None]
    selbond1            = (string)             [default: None]
    selbond2            = (string)             [default: None]
    targetcoordinate    = [float,float,float]  [default: None]
    remove_self         = (True/False)         [default: False]

    [Solvent]
    include_solvent     = (True/False)         [default: False]
    solvent_cutoff      = (float)              [default: None]
    solvent_selection   = (string)             [default: None]

    [Time]
    begintime           = (integer)            [default: 0]
    endtime             = (integer)            [default: None]
    dt                  = (integer)            [default: 1]
    skip                = (integer)            [default: 1]


A complete explanation of each option in the configuration file is available via the command:

    raijin -h

***Raijin*** has 3 calculations MODES:

* In ***ATOM*** mode, the coordinate of one atom will be tracked throughout the trajectory to serve as target point.
If more than 1 atom is provided in the selection, the center of geometry (COG) is used as target position.

* In ***BOND*** mode, the midpoint between 2 atoms will be tracked throughout the trajectory to serve as target
point. In ***BOND*** mode, the bond axis is used to calculate electric field alignment. By default, the bond axis is
define as **selbond1 ---> selbond2**.

* In ***COORDINATE*** mode, a list of [X,Y,Z] coordinates will serve as target point in all trajectory frames.


***IMPORTANT***:
* All selections must be compatible with MDAnalysis.
* If COORDINATE mode, make sure you have fixed translations and orientations in your trajectory.


## Raijin PyMOL plugin Usage

To install the PyMOL plugin, open PyMOL > Plugin Manager and click on "Install New Plugin" tab.
Load the ***Raijin*** plugin and use it via command-line within PyMOL. To usage instructions, read our FAQ.

Our plugin has 3 functions:

* **draw_bond_axis**: allows user to draw a vector illustrating the bond axis. By default: selection1 ---> selection2:
    draw_bond_axis resid 160 and name OG, resname LIG and name C1, color=yellow, gap=0.5

* **efield_point**: allows user to draw a vector illustrating the Electric field at a given coordinate point.
    efield_point [0,0,0], [-1.68931390e+02,1.64390616e+02,6.58994183e+00], color=magenta, scale=0.01

* **efield_bond**: allows user to draw a vector illustrating the Electric field at the midpoint between 2 atoms.
    efield_point resid 160 and name OG, resname LIG and name C1, [-1.68931390e+02,1.64390616e+02,6.58994183e+00], color=magenta, scale=0.01
