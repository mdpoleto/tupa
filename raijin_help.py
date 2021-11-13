#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#Marcelo D. Poleto
#APRIL 2022


header = """
	###########################################################################
	#                      Raijin - Electric Field calculator                 #
	#                            Raijin 1.0.0 (2022)                          #
	#                                                                         #
	#                    Marcelo D Poleto, Justin A Lemkul                    #
	#                                                                         #
	#                    Paper metadata here..............                    #
	#                                                                         #
	#            In case of bugs or suggestions, please contact us:           #
	#                              mdpoleto@vt.edu                            #
	#                                                                         #
	# LGPL-3.0                                                                #
	###########################################################################\n\n"""

template_content = """[Elecfield Selection]
# The atoms from which we calculate the electric field
sele_elecfield      = segid PROA

[Target Selection]
# Provide the target selection for the MODE of you choice
# e.g. if bond is used, then modebond1 and modebond2 must be defined.
mode                = ATOM    # ATOM or BOND or COORDINATE
selatom             = segid PROA and (resid 160 and name OG)
selbond1            = segid PROA and (resid 160 and name OG)
selbond2            = segid LIG and name C1
targetcoordinate    = [0,0,0]
remove_self         = True    # If target selection in sele_elecfield, remove it
						  # from sele_elecfield (True of False)

[Solvent]
include_solvent     = True    # or False
solvent_cutoff      = 20      # in Angstrom
solvent_selection   = segid TIP3

[Time]
begintime           = None    # begintime and endtime allow the user to evaluate just a part of the trajectory
endtime             = None
dt                  = 10      # Frequency of frames written in your trajectory (in picosecond)
skip                = 1       # How many frames do you want to skip?

# IMPORTANT:
# 1- All selections must be compatible with MDAnalysis
# 2- remove_self can be implicitely used by removing whatever contribution you
# do not want from sele_elecfield. Be smart about your selections.
# 3- the ATOM mode uses 1 atom to track its position througout trajectory and
# calculates the Efield at its position. Useful if you need the Efield sensed
# by a moving atom.
# 4- The BOND mode uses 2 atoms to calculate the Efield at a position equidistant
# to the 2 atoms selected. It uses the bond axis (default is selbond1 --> selbond2)
# to check Efield alignment, which is useful for catalysis engineering.
# 5- If COORDINATE mode, make sure you have fixed translations and
# orientations in your trajectory.
"""

# I NEED TO WORK ON THE HELP HERE!
help = """

  # Configuration File Inputs
  ###########################################################################
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



  # Configuration File Help
  ###########################################################################
  [Elecfield Selection]
  sele_elecfield      = selection of atoms that exert the electric field in the
                        calculation. Selection must be compatible with MDanalysis.

  [Target Selection]
  mode                = defines the tipe of calculation to be done. Possible
                        values are: ATOM, BOND and COORDINATE.

                        In ATOM mode, the coordinate of one atom will be tracked
                        throughout the trajectory to serve as target point.

                        In BOND mode, the midpoint between 2 atoms will be
                        tracked througout the trajectory to serve as target
                        point. In BOND mode, the bond axis is used to calculate
                        electric field alignment. By default, the bond axis is
                        define as selbond1 ---> selbond2.

                        In COORDINATE mode, a list of [X,Y,Z] coordinates will
                        serve as target point in all trajectory frames.

  selatom             = selection used in ATOM mode (compatible with MDanalysis)
  selbond1            = 1st selection used in BOND mode (compatible with MDanalysis)
  selbond2            = 2nd selection used in BOND mode (compatible with MDanalysis)
  targetcoordinate    = list of [X,Y,Z] coordinates (in Angstroms)
  remove_self         = True or False. Remove the contribution of atoms defined
                        in selatom or selbond1/selbond2 to the calculated electric
                        field. In practice, the behavior of remove_self = True is
                        identical of removing selatom or selbond1/selbond2 when
                        defining sele_elecfield.

  [Solvent]
  include_solvent     = True or False. Either account for the contribution of
                        solvent atoms to the exerted electric field or not.
  solvent_cutoff      = cutoff radius used to select solvent atoms in relation
                        to the selection used in sele_elecfield.
  solvent_selection   = selection of solvent entities to be accounted for in the
                        electric field calculation.

  [Time]
  begintime           = begintime and endtime (in ps) allow the user to analyze
                        only a chunk of the entire trajectory.
  endtime             =
  dt                  = frequency in which frames were written in your trajectory.
                        dt is used to convert frame number into simulation time.
  skip                = allow user to skip an amount of frames when analyzing
                        the trajectory.
"""
