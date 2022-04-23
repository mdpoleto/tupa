[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)
[![Twitter Follow](https://img.shields.io/twitter/follow/mdpoleto?style=social)](https://twitter.com/mdpoleto)

# **TUPÃ**: Electric field analyses for molecular simulations

<img src="https://github.com/mdpoleto/tupa/blob/main/Figures/TUPÃ_LOGO.png">

## What is TUPÃ?
**TUPÃ** (pronounced as [*tu-pan*](https://translate.google.com/?hl=pt-BR&sl=pt&tl=en&text=tup%C3%A3&op=translate)) is a python algorithm that employs MDAnalysis engine to calculate electric fields at any point inside
the simulation box throughout MD trajectories. **TUPÃ** also includes a PyMOL plugin to visualize electric
field vectors together with molecules.

Required packages:

* MDAnalysis >= 1.0.0
* Python     >= 3.x
* Numpy      >= 1.2.x

------------------------------
## Installation instructions

First, make sure you have all required packages installed. For MDAnalysis installation procedures, [click here](https://www.mdanalysis.org/pages/installation_quick_start/).

After, just clone this repository into a folder of your choice:
```
git clone https://github.com/mdpoleto/tupa.git
```

To easily use **TUPÃ**, export the pathway to the executable *TUPA.py* in your ~/.bashrc:
```
export PATH=/path/to/the/cloned/repository/:$PATH
```


## TUPÃ Usage
**TUPÃ** calculations are based on parameters that are provided via a configuration file,
which can be obtained via the command:
```
TUPA.py -template config.conf
```

The configuration file usually contains:
```
[Environment Selection]
sele_environment    = (string)             [default: None]

[Probe Selection]
mode                = (string)             [default: None]
selatom             = (string)             [default: None]
selbond1            = (string)             [default: None]
selbond2            = (string)             [default: None]
probecoordinate     = [float,float,float]  [default: None]
file_of_coordinates = (pathway to file)    [default: None]
remove_self         = (True/False)         [default: False]
remove_cutoff       = (float)              [default: 1 A ]

[Solvent]
include_solvent     = (True/False)         [default: False]
solvent_cutoff      = (float)              [default: 10 A]
solvent_selection   = (string)             [default: None]

[Time]
dt                  = (integer)            [default: 1]

[Box Info]
redefine_box        = Whether or not provide explicit box dimension information.
boxdimensions       = Box dimension information [A,B,C,Alpha,Beta,Gamma]. A,B
                      and C are the edge lengths (in Angstrom). Alpha, Beta
                      and Gamma are the box internal angles (in degrees)
```

A complete explanation of each option in the configuration file is available via the command:
```
TUPA.py -h
```

**TUPÃ** has 4 calculations MODES:

* In ``ATOM`` mode, the coordinate of one atom will be tracked throughout the trajectory to serve as probe point.
If more than 1 atom is provided in the selection, the center of geometry (COG) is used as probe position. An example
is provided [HERE](https://github.com/mdpoleto/tupa/tree/main/Examples/ATOM).

* In ``BOND`` mode, the midpoint between 2 atoms will be tracked throughout the trajectory to serve as probe
point. In this mode, the bond axis is used to calculate electric field alignment. By default, the bond axis is
define as ```selbond1 ---> selbond2```. An example is provided [HERE](https://github.com/mdpoleto/tupa/tree/main/Examples/BOND).

* In ``COORDINATE`` mode, a [X,Y,Z] coordinate will serve as probe point in all trajectory frames.
An example is provided [HERE](https://github.com/mdpoleto/tupa/tree/main/Examples/COORDINATE).

* In ``LIST`` mode, a list of [X,Y,Z] coordinates will serve as probe points, one for each trajectory frame.

**IMPORTANT**:
* All selections must be compatible with MDAnalysis syntax.
* TUPÃ was designed to work with ```ORTHORHOMBIC``` box types. We are working to support for rhombic dodecahedron and truncated octahedron boxes.
* Trajectories MUST be re-imaged before running TUPÃ. *Make sure your probe is well centered in the box*.
* Molecules in ```solvent_selection``` beyond the PBC are re-imaged. This is achieved by applying the ```around``` selection feature in MDAnalysis and properly shifting the coordinates.
* If using COORDINATE mode, be mindful that our code does not account for rotations and translations of coordinates. Be mindful of the coordinate selection.


## TUPÃ PyMOL Plugin (pyTUPÃmol)

<img src="https://github.com/mdpoleto/tupa/blob/main/Figures/pyTUPÃmol_example.png" width="680">

**pyTUPÃmol** is a PyMOL plugin to plot electric field vectors alongside other molecules. By definition, an arrow CGO object is created starting at a given coordinate [X,Y,Z] and has the size of the magnitude of the provided electric field.

Our plugin has 3 functions that can be called via command line within PyMOL:

* **efield_point**: create a vector at a given atom or set of coordinates.
```
efield_point segid LIG and name O1, efield=[-117.9143, 150.3252, 86.5553], scale=0.01, color="red", name="efield_OG"
```

* **efield_bond**: create a vector midway between 2 selected atoms.
```
efield_point resname LIG and name O1, resname LIG and name C1, efield=[-94.2675, -9.6722, 58.2067], scale=0.01, color="blue", name="efield_OG-C1"
```

* **draw_bond_axis**: create a vector representing the axis between 2 atoms.
```
draw_bond_axis resname LIG and name O1, resname LIG and name C1, gap=0.5, color="gray60", name="axis_OG-C1"
```

To install **pyTUPÃ** plugin in PyMOL, click on Plugin > Plugin Manager and then "Install New Plugin" tab. Choose the ```pyTUPÃ.py``` file and click Install.

--------------------------
## Citing TUPÃ

If you use **TUPÃ** in a scientific publication, we would appreciate citations to the following paper:

Marcelo D. Polêto, Justin A. Lemkul. *TUPÃ: Electric field analysis for molecular simulations*, 2022.

Bibtex entry:
```
@article{TUPÃ2022,
    author = {Pol\^{e}to, M D and Lemkul, J A},
    title = {TUPÃ: Electric field analyses for molecular simulations},
    journal = {Journal of Computational Chemistry},
    volume = {n/a},
    number = {n/a},
    pages = {},
    keywords = {electric field, electrostatics, force fields, molecular dynamics, molecular mechanics},
    doi = {https://doi.org/10.1002/jcc.26873},
    url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.26873},
    eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/jcc.26873},
    abstract = {Abstract We introduce TUPÃ, a Python-based algorithm to calculate and analyze electric fields in molecular simulations. To demonstrate the features in TUPÃ, we present three test cases in which the orientation and magnitude of the electric field exerted by biomolecules help explain biological phenomena or observed kinetics. As part of TUPÃ, we also provide a PyMOL plugin to help researchers visualize how electric fields are organized within the simulation system. The code is freely available and can be obtained at https://mdpoleto.github.io/tupa/.}
}
```
## Why TUPÃ?

In Brazilian folklore, Tupã is considered a "manifestation of God in the form of thunder". To know more, refer to [this](https://en.wikipedia.org/wiki/Tup%C3%A3_(mythology)).

## Contact information
E-mail: mdpoleto@vt.edu / jalemkul@vt.edu

<a href="https://trackgit.com">
<img src="https://us-central1-trackgit-analytics.cloudfunctions.net/token/ping/l02lg00zonn9v19irctl" alt="trackgit-views" />
</a>
