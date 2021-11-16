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


## Installation instructions
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
