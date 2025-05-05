#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#Marcelo D. Poleto
#Apr 2023


header = """
	###########################################################################
	#                 TUPÃ - Electric Field analyses algorithm                #
	#                            TUPÃ v1.6.0 (2022)                           #
	#                                                                         #
	#                    Marcelo D Poleto, Justin A Lemkul                    #
	#                                                                         #
	#          TUPÃ: Electric field analyses for molecular simulations        #
	#                  Journal of Computational Chemistry (2022)              #
	#                   url: https://doi.org/10.1002/jcc.26873                #
	#                                                                         #
	#            In case of bugs or suggestions, please contact us:           #
	#                              mdpoleto@vt.edu                            #
	#                                                                         #
	# LGPL-3.0                                                                #
	###########################################################################\n\n"""

template_content = """[Environment Selection]
# The atoms from which we calculate the electric field
sele_environment    = resid 1:200


[Probe Selection]
# Provide the probe selection for the MODE of you choice
# e.g. if bond is used, then selbond1 and selbond2 must be defined.
mode                = ATOM    # ATOM or BOND or COORDINATE or LIST
selatom             = resname LIG and name C5
selbond1            = resname LIG and name C5
selbond2            = resname LIG and name C1
probecoordinate     = [0,0,0]
file_of_coordinates = /path/to/file/of/list/of/probe/coordinates

remove_self         = True    # For COORDINATE and LIST modes only, whether remove the
                              # contribution of self within a cutoff of the coordinate
remove_cutoff       = 1       # in Angstrom


[Solvent]
include_solvent     = True    # or False
solvent_cutoff      = 15      # in Angstrom
solvent_selection   = resname TIP3 or resname CL

[Time]
dt                  = 10      # Frequency of frames written in your trajectory (in picosecond)

[Box Info]
redefine_box        = False   # Your box SHOULD have dimensions for each frame!
                              # This is a workaround in case you don't.
boxdimensions       = [float,float,float,float,float,float]

# IMPORTANT:
# 1- All selections must be compatible with MDAnalysis
# 2- remove_self only works in COORDINATE and LIST modes. For ATOM or BOND mode, analogous
# behavior can be created by removing whatever contribution you do not want from
# sele_environment. Be smart about your selections.
# 3- the ATOM mode uses 1 atom to track its position througout trajectory and
# calculates the Efield at its position. If more than 1 atom is provided in selection,
# the center of geometry (COG) is used as target position.
# 4- The BOND mode uses 2 atoms to calculate the Efield at a position equidistant
# to the 2 atoms selected. It uses the bond axis (default is selbond1 --> selbond2)
# to check Efield alignment, which is useful for catalysis engineering. Note that
# a bond dipole direction is defined to be from positive to negative.
# 5- If COORDINATE mode, make sure you have fixed translations and
# orientations in your trajectory.
# 6- Your box should have dimensions for each frame. If you these are missing, we
# suggest to go back to the step where you reoriented your trajectory and double
# check your steps there.
"""

help = """
  # Configuration File Inputs
  ###########################################################################
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
  redefine_box        = (True/False)         [default: False]
  boxdimensions       = [float,float,float,float,float,float] [default: None]


  # Configuration File Help
  ###########################################################################
  [Environment Selection]
  sele_environment    = selection of atoms that exert the electric field in the
                        calculation. Selection must be compatible with MDanalysis.

  [Probe Selection]
  mode                = defines the tipe of calculation to be done. Possible
                        values are: ATOM, BOND, COORDINATE and LIST.

                        In ATOM mode, the coordinate of one atom will be tracked
                        throughout the trajectory to serve as probe point. If
                        more than 1 atom is provided in the selection, the
                        center of geometry (COG) is used as probe position.

                        In BOND mode, the midpoint between 2 atoms will be
                        tracked throughout the trajectory to serve as probe
                        point. In BOND mode, the bond axis is used to calculate
                        electric field alignment. By default, the bond axis is
                        define as selbond1 ---> selbond2.

                        In COORDINATE mode, a [X,Y,Z] coordinate will
                        serve as probe point in all trajectory frames.

                        In LIST mode, a list of [X,Y,Z] coordinates will
                        serve as probe points, one for each trajectory frame.

  selatom             = selection used in ATOM mode (compatible with MDanalysis)
  selbond1            = 1st selection used in BOND mode (compatible with MDanalysis)
  selbond2            = 2nd selection used in BOND mode (compatible with MDanalysis)
  probecoordinate     = [X,Y,Z] coordinate (in Angstroms)
  file_of_coordinates = pathway to file of a list of coordinates to be tracked as probes
  remove_self         = True or False. Only works for COORDINATE mode. It removes
                        the contribution of atoms defined in sele_environment within
                        within a determined cutoff to the electric field.
  remove_cutoff       = cutoff used to determine wether to remove the self
                        contribution or not.


  [Solvent]
  include_solvent     = True or False. Either account for the contribution of
                        solvent atoms to the exerted electric field or not.
  solvent_cutoff      = cutoff radius used to select solvent atoms in relation
                        to the selection used in sele_environment.
  solvent_selection   = selection of solvent entities to be accounted for in the
                        electric field calculation.

  [Time]
  dt                  = frequency in which frames were written in your trajectory.
                        dt is used to convert frame number into simulation time.

  [Box Info]
  redefine_box        = Whether or not to overwrite box dimension information.
  boxdimensions       = Box dimension information [A,B,C,Alpha,Beta,Gamma]. A,B
                        and C are the edge lengths (in Angstrom). Alpha, Beta
                        and Gamma are the box internal angles (in degrees)
"""
