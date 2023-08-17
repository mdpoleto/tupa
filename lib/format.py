#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#Marcelo D. Poleto
#Apr 2023

import os
from lib.utils import *


##############################################################
# Opening output files
#
# - run_info.log
# - ElecField.dat
# - ElecField_per_residue.dat
# - Spatial_Deviation.dat
# - probe_info.dat
#
# if bond mode:
# - ElecField_proj_onto_bond.dat
# - ElecField_alignment.dat
# - bond_axis_info.dat
#
##############################################################

def print_run_info(parsed_configuration, top_file, traj_file, outdir):

	##############################
	# Defining variables
	sele_elecfield    = parsed_configuration.sele_elecfield

	mode              = parsed_configuration.mode
	if mode == "atom":
	    selatom = parsed_configuration.selatom
	    remove_cutoff = 0
	elif mode == "bond":
	    selbond1 = parsed_configuration.selbond1
	    selbond2 = parsed_configuration.selbond2
	    remove_cutoff = 0
	elif mode == "coordinate":
	    probecoordinate = parsed_configuration.probecoordinate
	    remove_self = parsed_configuration.remove_self
	    remove_cutoff = parsed_configuration.remove_cutoff
	elif mode == "list":
	    file_of_coordinates = parsed_configuration.file_of_coordinates
	    remove_self = parsed_configuration.remove_self
	    remove_cutoff = parsed_configuration.remove_cutoff

	include_solvent   = parsed_configuration.include_solvent
	if include_solvent:
	    solvent_selection = parsed_configuration.solvent_selection
	    solvent_cutoff = parsed_configuration.solvent_cutoff

	dt = parsed_configuration.dt

	redefine_box = parsed_configuration.redefine_box
	if redefine_box:
	    boxdimensions = parsed_configuration.boxdimensions


	##############################
	# Defining probe string for each calculation mode
	probestr = ""
	if mode == "atom":
		str1 = 'selatom             = {}'.format(selatom)
		probestr += str1
	elif mode == "bond":
		str1 = 'selbond1            = {}'.format(selbond1)
		str2 = 'selbond2            = {}'.format(selbond2)
		probestr += str1 +'\n'+str2
	elif mode == "coordinate":
		str1 = 'probecoordinate     = {}'.format(probecoordinate)
		probestr += str1
	elif mode == "list":
		str1 = 'file_of_coordinates = {}'.format(file_of_coordinates)
		probestr += str1+'\n'
		if remove_self == True:
			str2 = 'remove_self         = {}'.format(remove_self)
			str3 = 'remove_cutoff       = {}'.format(remove_cutoff)
			probestr += str2 +'\n' + str3
	solvstr = ""
	if include_solvent == True:
		str1 = 'solvent_selection   = {}'.format(solvent_selection)
		str2 = 'solvent_cutoff      = {}'.format(solvent_cutoff)
		solvstr += str1 +"\n"+str2
	boxstr = ""
	if redefine_box == True:
		str1 = 'boxdimensions       = {}'.format(boxdimensions)
		boxstr += str1


	##############################################################
	# Write the parameters to a file

	if not os.path.exists(outdir):
		os.system("mkdir " + outdir)

	# Writing run info into a file
	out = open(os.path.join(outdir,"run_info.log"),'w')

	runinfo = """########################################################
>>> Parameters used to run TUPÃ:

Topology file       = {top}
Trajectory file     = {traj}

[Environment Selection]
sele_environment    = {environment}

[Probe Selection]
mode                = {mode}
{probestr}

[Solvent]
include_solvent     = {include_solv}
{solvstr}

[Time]
dt                  = {dt}

[Box Info]
redefine_box        = {box}
{boxstr}
########################################################
	""".format(top = top_file, traj = traj_file,
	           environment = sele_elecfield,
	           mode=mode,
	           probestr=probestr,
	           include_solv=include_solvent,
	           solvstr=solvstr,
	           dt=dt,
	           box=redefine_box,
	           boxstr=boxstr)

	out.write(runinfo)
	out.close()

	print(runinfo)


def open_outputs(outdir):
	out = open(os.path.join(outdir,"ElecField.dat"),'w')
	out.write("""@    title "Electric Field"\n@    xaxis  label "Time (ps)"\n@    yaxis  label "MV/cm"\n""")
	out.write("#time             Magnitude      Efield_X       Efield_Y       Efield_Z\n")
	out.write("@type xy\n")

	outres = open(os.path.join(outdir,"ElecField_per_residue.dat"),'w')
	outres.write("""@    title "Magnitude of Electric Field"\n@    xaxis  label "Residues"\n@    yaxis  label "MV/cm"\n""")
	outres.write("#time              Ef_per_res     Std.Dev        Alignment%     Alignment_Std\n")
	outres.write("@type xydy\n")

	outangle = open(os.path.join(outdir,"Spatial_Deviation.dat"), "w")
	outangle.write("#time   Angle(Efield(t), avg_field)   Projection(Efield(t), avg_field)   Alignment(Efield(t), avg_field)\n")

	outprobe = open(os.path.join(outdir,"probe_info.dat"), "w")
	outprobe.write("#time     probe_coordinates\n")

	return out, outres, outangle, outprobe


def open_outputs_bond(outdir):
	outproj = open(os.path.join(outdir,"ElecField_proj_onto_bond.dat"),'w')
	outproj.write("""@    title "Efield_Projection"\n@    xaxis  label "Time (ps)"\n@    yaxis  label "MV/cm"\n""")
	outproj.write("#time             Magnitude       Efield_X       Efield_Y      Efield_Z        Direction\n")
	outproj.write("@type xy\n")

	outalig = open(os.path.join(outdir,"ElecField_alignment.dat"),'w')
	outalig.write("""@    title "Efield_Projection"\n@    xaxis  label "Time (ps)"\n@    yaxis  label "Alignment_percentage"\n""")
	outalig.write("#time             Angle           Alignment%\n")
	outalig.write("@type xy\n")

	outresbond = open(os.path.join(outdir,"Residue_alignment_to_bond.dat"),'w')
	outresbond.write("""@    title "Magnitude of Electric Field"\n@    xaxis  label "Residues"\n@    yaxis  label "MV/cm"\n""")
	outresbond.write("#time              Efproj_res     Std.Dev        AVG_align%     Std_align%\n")
	outresbond.write("@type xydy\n")

	outrbond = open(os.path.join(outdir,"bond_axis_info.dat"), "w")
	outrbond.write("#time       rbondvec_in_meters                                   rbond_magnitude             rbond_unit_vector\n")

	return outproj, outalig, outrbond, outresbond


##############################################################
# Formatting output lines
#
# - ElecField.dat
# - ElecField_per_residue.dat
# - Spatial_Deviation.dat
#
# if bond mode:
# - ElecField_proj_onto_bond.dat
# - ElecField_alignment.dat
# - bond_axis_info.dat
#
##############################################################
def fmt_efield_out_line(field_array, t=None, lastline = False):
	if not lastline:
		string = "{a:<12}{b: 15.6f}{c: 15.6f}{d: 15.6f}{e: 15.6f}\n".format(a=t, b=field_array[0], c=field_array[1], d=field_array[2],e=field_array[3])
	else:
		string = "{a: 15.6f}{b: 15.6f}{c: 15.6f}{d: 15.6f}\n".format(a=field_array[0], b=field_array[1], c=field_array[2], d=field_array[3])
	return string


def fmt_res_out_line(field_array):
	string = "{a:<12}{b: 15.6f}{c: 15.6f}{d: 15.6f}{e: 15.6f}\n".format(a=field_array[0], b=field_array[1], c=field_array[2], d=field_array[3],e=field_array[4])
	return string


def fmt_spatialdev_out_line(field_array,t=None, lastline=False):
	if lastline == False:
		string = "{a:<12}{b: 15.6f}{c: 30.6f}{d: 33.6f}\n".format(a=t, b=field_array[0], c=field_array[1], d=field_array[2])
	else:
		string = "{a: 15.6f}{b: 30.6f}{c: 33.6f}\n".format(a=field_array[0], b=field_array[1], c=field_array[2])
	return string


def fmt_proj_out_line(field_array, t=None, lastline = False, optarr=None):
	if not lastline:
		string = "{a:<12}{b: 15.6f}{c: 15.6f}{d: 15.6f}{e: 15.6f}{f: 15.6f}\n".format(a=t, b=field_array[0], c=field_array[1], d=field_array[2], e=field_array[3], f=optarr[0])
	else:
		string = "{a: 15.6f}{b: 15.6f}{c: 15.6f}{d: 15.6f}\n".format(a=field_array[0], b=field_array[1], c=field_array[2], d=field_array[3])
	return string


def fmt_align_out_line(field_array, lastline = False):
	if not lastline:
		string = "{a:<12}{b: 15.6f}{c: 15.6f}\n".format(a=field_array[0], b=field_array[1], c=field_array[2])
	else:
		string = "{a: 15.6f}{b: 15.6f}\n".format(a=field_array[0], b=field_array[1])
	return string


def fmt_rbond_info(field_array):
	rvx, rvy, rvz = field_array[1]
	rhx, rhy, rhz = field_array[2]
	list1 = "[" + str(rvx) + "," + str(rvy) + "," + str(rvz) + "]"
	list2 = "[" + str(rhx) + "," + str(rhy) + "," + str(rhz) + "]"
	string = "{a:<12}{b:50s}{c: 11.6f}                    {d:<50s}\n".format(a=field_array[0], b=list1, c=field_array[3], d=list2)
	return string

def fmt_probe_position(time, array):
	string = "[" + str(array[0]) + ", " + str(array[1]) + ", " + str(array[2]) + "]"
	string = str(time).ljust(10,' ') +  string.ljust(60,' ') + "\n"
	return string
