#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#Marcelo D. Poleto
#MAR 2023

import os
import numpy as np


##############################################################
# Defining analytical functions
#
# - mag()
# - unit_vector()
# - alignment()
# - projection()
# - angle_between()
# - calc_ElectricField()
# - pack_around()
#
##############################################################

def mag(vector):
    """ Returns the magnitude of the vector.  """
    mag = np.linalg.norm(vector)
    return mag


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def alignment(value1,value2):
    """ Returns % alignment between 2 values """
    # V2 is the total reference
    aligned = abs(value1/value2)*100
    return aligned


def projection(v1,v2):
    """ Returns projection of vector v1 onto vector v2 """
    #projection_u_on_v = (np.dot(u, v)/np.dot(v, v))*v
    proj = v2*(np.dot(v1, v2)/np.dot(v2,v2))
    return proj


def angle_between(v1, v2):
    """ Returns angle (in degrees) between vector v1 and vector v2 """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    angle_deg = np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))
    return angle_deg


def calc_ElectricField(atom,refposition):
    """Calculate electric field exerted by 'atom' at 'refposition'
    Parameters
    ----------
    atom: Atom MDA object
        input DataFrame
    refposition: np.array 
        numpy array of shape (N,3) with X Y Z coordinates
    Returns
    -------
    Ef: np.array (N,3)
        numpy array of shape (N,3) with X Y Z electric field vectors
    """
   
    #Epsilon = 8.8541878128e-12 C**2/N*(m**2)
    k =  8987551792.261173  # in (N*m**2)/C**2, pre-computed from  1 / (4 * np.pi * Epsilon) 

    rvec = (refposition * 1e-10) - (atom.position * 1e-10)  # Converting A to meters

    rmag = mag(rvec)
    rhat = unit_vector(rvec)

    # convert charge (elementary charge unit to Coulomb)
    # 1 elementary charge unit = 1.60217733*(10**(-19))
    charge = np.round(atom.charge, 6) * 1.60217733e-19 # convert to Coulomb

    # Calculate Electric Field in N/C
    # Ef = k*Q/r2; charge = coulomb ; r = meter
    Ef = rhat *( k * charge / rmag ** 2)
    Ef *= 1e-8 #Convert to MegaVolt per centimeter (usual unit in literature)
    # 1 MV/cm = 10**8 V/m (the SI unit), 1.944*10**(−4) atomic units
    # 10 MV/cm = 1 V/nm = 0.1 V/Angstrom

    return Ef


def pack_around(universe, atom_group, center):
    """ Transpose coordinates in a rectangular PBC box """
    tmp_group = atom_group # we will update tmp_group below
    atom_group.pack_into_box()
    positions = atom_group.positions.copy() # working with a tmp copy of positions
    sub = positions - center

    # Get the box for the current frame
    box = universe.dimensions[:3]
    culprits = np.where(np.sqrt(sub**2) > box / 2)
    # Actually translate the coordinates.
    positions[culprits] -= (universe.dimensions[culprits[1]]
                            * np.sign(sub[culprits]))
    tmp_group.positions = positions

    return tmp_group


def check_box_vectors(universe, config):
    """ Check whether trajectory file has box dimension information and
        redefine it (if it has been requested in the configuration file) """

    if universe.dimensions is not None:
        boxangles = universe.dimensions[3:]
        checkangles = [i for i in boxangles if float(i) != 90.]
        if len(checkangles) == 0:
            if universe.dimensions[0] == 1 and universe.dimensions[1] == 1 and universe.dimensions[2] == 1:
                if config.redefine_box == True:
                    print("\n>>> Redefining box dimensions to:", config.boxdimensions)
                    always_redefine_box_flag = True
                else:
                    raise ValueError(""">>> ERROR: Your trajectory does not contain information regarding box size. Provide them in the configuration file!\n""")
            elif universe.dimensions[0] == 0 and universe.dimensions[1] == 0 and universe.dimensions[2] == 0:
                print("\n>>> Redefining box dimensions to: ", config.boxdimensions)
                always_redefine_box_flag = True
            else:
                if config.redefine_box == True:
                    print("\n>>> Redefining box dimensions to:", config.boxdimensions)
                    always_redefine_box_flag = True
                else:
                    always_redefine_box_flag = False
        else:
            raise ValueError(""">>> WARNING: Your box does not seem to be orthorhombic. We are currently working to support non-rectangular boxes.\n""")
            always_redefine_box_flag = False
    elif universe.dimensions is None and config.redefine_box is True:
        always_redefine_box_flag = True
    else:
        raise ValueError(""">>> ERROR: Your trajectory does not contain information regarding box size. Provide them in the configuration file!\n""")

    return always_redefine_box_flag


def dump_coor_time(universe, time, dumptime, enviroment_selection, outdir):
    dumptime = np.array(dumptime, dtype='float64')
    if float(time) in dumptime:
        framefile = os.path.join(outdir,"frame_" + str(time) + "ps.pdb")
        envfile   = os.path.join(outdir,"environment_" + str(time) + "ps.pdb")
        print("   >>> Dumping frame (Time = " + str(time) + " ps)! Check " + framefile)
        dump_sel = universe.select_atoms("all")
        dump_sel.write(framefile)
        dump_sel2 = enviroment_selection.select_atoms("all")
        dump_sel2.write(envfile)


def create_probe_selection(universe, config):
    """Creating the probe selection. Being verbose so users can better
       understand what is being defined/done. """

    #####################################
    # Verbose output for solvent inclusion in calculation. Selection is done within the trajectory loop
    if config.include_solvent == True:
        print(">>> Including solvent molecules within a radius of " + str(config.solvent_cutoff) + " Angstroms in the calculation...")
    else:
        print(">>> Excluding solvent molecules from the calculation...")


    if config.mode == "atom":
        try:
            probe_selection = universe.select_atoms(config.selatom, updating=True)
            if not probe_selection:
                raise ValueError(">>> ERROR: PROBE selection (from selatom) is empty. Check your selection!\n")
        except:
            raise ValueError(">>> ERROR: PROBE selection (from selatom) could not be parsed. Check your selection!\n")

        if len(probe_selection) > 1:
            print("\n>>> Running in ATOM mode!")
            print(">>> Probe atoms (using center of geometry) = "+ str(probe_selection.atoms) + "\n")
        else:
            print("\n>>> Running in ATOM mode!")
            print(">>> Probe atom = "+ str(probe_selection.atoms[0].name) +"-"+ str(probe_selection.resnames[0]) + str(probe_selection.resnums[0]) + "\n")

        return probe_selection

    elif config.mode == "bond":
        try:
            bond1 = universe.select_atoms(config.selbond1, updating=True)
            if not bond1:
                raise ValueError(">>> ERROR: PROBE selection (from selbond1) is empty. Check your selection!\n")
        except:
            raise ValueError(">>> ERROR: PROBE selection (from selbond1) could not be parsed. Check your selection!\n")
        try:
            bond2 = universe.select_atoms(config.selbond2, updating=True)
            if not bond2:
                raise ValueError(">>> ERROR: PROBE selection (from selbond2) is empty. Check your selection!\n")
        except:
            raise ValueError("\n>>> ERROR: PROBE selection (from selbond2) could not be parsed. Check your selection!\n")

        probe_selection = bond1 + bond2
        print("\n>>> Running in BOND mode!")
        print(">>> Probe atoms = " + str(probe_selection.atoms))
        print(">>> Bond axis direction = " + str(bond1.atoms[0].name) +"-"+ str(bond1.resnames[0]) + str(bond1.resnums[0]) + " --> " + str(bond2.atoms[0].name) +"-"+ str(bond2.resnames[0]) + str(bond2.resnums[0]) + "\n")

        return probe_selection

    elif config.mode == "coordinate":
        probe_selection = config.probecoordinate
        print("\n>>> Running in COORDINATE mode!")
        print(">>> Probe XYZ position = " +  str(probe_selection))
        if config.remove_self == True:
            print(">>> Removing self contribution within a radial cutoff of " +  str(config.remove_cutoff) + " Angstroms\n")

        return probe_selection

    elif config.mode == "list":
        print("\n>>> Running in LIST mode!")
        # If mode == LIST, parse the file of coordinates
        loc = open(config.file_of_coordinates, 'r')
        probe_selection = []
        while True:
            coorline = loc.readline()
            if not coorline: break
            if coorline[0] != "#" and coorline[0] != "@":
                if len(coorline.split(",")) == 3:
                    try:
                        listcoorX = float(coorline.split(",")[0])
                        listcoorY = float(coorline.split(",")[1])
                        listcoorZ = float(coorline.split(",")[2])
                    except:
                        raise ValueError(""">>> ERROR: File of coordinates could not be parsed! Expecting floats or integers (X, Y, Z). Check your inputs!\n""")
                    
                    tmpcoorlist = np.array([listcoorX, listcoorY, listcoorZ])
                    probe_selection.append(tmpcoorlist)
                else:
                    raise ValueError(""">>> ERROR: File of coordinates could not be parsed! Expecting 3 columns (X, Y, Z). Check your inputs!\n""")

        # Sanity check: list of coordinates should have the same length as the trajectory
        if len(probe_selection) != len(universe.trajectory):
            raise ValueError(""">>> ERROR: File of coordinates and trajectory have different lengths (""" + str(len(probe_selection)) + """ lines and """ + str(len(universe.trajectory)) + """ frames)!\n""")

        print(">>> Probe XYZ positions = " + str(len(probe_selection)) + """ XYZ coordinates""")
        if config.remove_self == True:
            print(">>> Removing self contribution within a radial cutoff of " +  str(config.remove_cutoff) + " Angstroms\n")

        return probe_selection

    else:
        raise ValueError(">>> MODE should be 'bond', 'atom', 'coordinate' or 'list'!")


def update_probe_position(universe, config, probe_selection):
    """Updating the probe position. No need to be verbose """
    if config.mode == "atom":
        if len(probe_selection) > 1:
            refposition = probe_selection.center_of_geometry()
        else:
            refposition = probe_selection.atoms[0].position

        return refposition

    elif config.mode == "bond":
        bond1 = universe.select_atoms(config.selbond1, updating=True)
        bond2 = universe.select_atoms(config.selbond2, updating=True)
        refposition = (bond1.atoms[0].position + bond2.atoms[0].position)/2 # midway for both atoms

        return refposition

    elif config.mode == "coordinate":
        refposition = probe_selection

        return refposition

    elif config.mode == "list":
        refposition = probe_selection[universe.trajectory.ts.frame] # update refposition from list
        return refposition


def update_environment(universe, config, elecfield_selection, refposition):
    """ Update environment to include solvent selection around a radius """
    # We incorporate the solvent selection around the probe here for each mode
    if config.include_solvent == True:
        tmp_selection = universe.select_atoms("(byres point " + str(refposition[0]) + " " + str(refposition[1]) + " " + str(refposition[2]) + " " + str(config.solvent_cutoff) + ") and (" + config.solvent_selection + ")", periodic=True)

        # translate molecules within selection that are beyond the PBC and incorporate into the environment
        tmp_selection = pack_around(universe, tmp_selection, refposition)
        enviroment_selection = elecfield_selection + tmp_selection
    else:
        enviroment_selection = elecfield_selection


    # Evaluate self_contribution removal
    # selects all atoms from environment within a cutoff of a point in space
    # (point X Y Z cutoff)
    # Such approach is only available to COORDINATE and LIST modes.
    # ATOM and BOND modes can achieve the same behavior by adjusting the
    # environment selection accordingly.
    self_contribution = enviroment_selection.select_atoms("byres point  " + str(refposition[0]) + " " + str(refposition[1]) + " " + str(refposition[2]) + " " + str(config.remove_cutoff) + " ", periodic=True)

    if len(self_contribution.atoms) > 0:
        if config.remove_self == True:
            print(""">>> WARNING! Removing self contribution of: """ + str(self_contribution.atoms))
            # Remove self_contribution for this frame
            enviroment_selection = enviroment_selection - self_contribution
        else:
            print(""">>> WARNING! Some atoms are closer than """ + str(config.remove_cutoff) + """ A : """ + str(self_contribution.atoms) + """ atoms""")

    return enviroment_selection


def calculate_Efprojection(universe, totalEf, config):
    """ Calculate Electric Field projection onto bond"""
    # Calculate efield projection
    bond1 = universe.select_atoms(config.selbond1, updating=True)
    bond2 = universe.select_atoms(config.selbond2, updating=True)
    rbond_vec = (bond2.atoms[0].position - bond1.atoms[0].position)* 1e-10 # convert from Angstrom to meter
    Efproj = projection(totalEf,rbond_vec)

    # Calculate projection direction (either 1 or -1)
    angle_deg = angle_between(totalEf,rbond_vec)
    proj_direction = np.cos(angle_deg)/abs(np.cos(angle_deg))

    return Efproj, proj_direction, angle_deg, rbond_vec


def add_to_dict(dict, key, array):
    """ Insert an array into a dictionary accounting for existing keys """
    if key in dict.keys():
        prev_array = dict[key]
        new_array = np.vstack((prev_array, array))
        dict[key] = new_array
    else:
        dict[key] = array

    return dict


class Results:
    def __init__(self):
        self.tmp_dict_res = {}
        self.res_contribution_per_frame = {}

        self.efield_timeseries = {}
        self.spatialdev = []
        self.resEFalignment_total = {}

        self.projection_timeseries = {}
        self.projalignment = {}
        self.resEFalignment_bond_per_frame = {}
