#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#Marcelo D. Poleto
#Apr 2023

import sys, os, argparse, timeit
import MDAnalysis as mda
import numpy as np
import configparser as cp
import warnings
sys.dont_write_bytecode = True
warnings.filterwarnings("ignore", message="Found no information for attr:")
warnings.filterwarnings("ignore", message="Found missing chainIDs")
np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

#####################################
# import methods and functions
from lib import tupa_help
from lib import format
from lib.utils import *
from lib.config_parser import Configuration
#####################################

#####################################
# import tools for parallelization
import multiprocessing
from multiprocessing import Pool
from functools import partial
#####################################


# Parse user input and options
ap = argparse.ArgumentParser(description=tupa_help.header,
                             formatter_class=argparse.RawDescriptionHelpFormatter,
                             usage="",
                             add_help=False)
ap._optionals.title = 'Options'
ap.add_argument('-h', '--help', action="store_true", help='Show this message and exit')
ap.add_argument('-top', type=str, default=None, required=False, metavar='',
                help='Topology file (.psf, .tpr, .prmtop, etc.)')
ap.add_argument('-traj', type=str, default=None, required=False, metavar='',
                help='Trajectory file (.dcd, .xtc, etc.)')
ap.add_argument('-outdir', type=str, default="Tupã_results/", required=False, metavar='',
                help='Output folder in which results will be written (default: Tupã_results/)')
ap.add_argument('-config', type=str, default=None, required=False, metavar='',
                help='Input configuration file (default: config.dat)')
ap.add_argument('-template', type=str, default=None, required=False, nargs="?", metavar='',
                help='Create a template for input configuration file')
ap.add_argument('-dumptime', type=int, default=None, required=False, metavar='', nargs='+',
                help='Choose times (in ps) to dump the coordinates from.')
ap.add_argument('-ncpus', type=int, default=1, required=False, metavar='',
                help='Number of CPUS to run calculations on (default: 1).')

cmd = ap.parse_args()
start = timeit.default_timer()

###############################################################################
# Printing help info and template
if cmd.help is True:
    ap.print_help()
    print(tupa_help.help)
    sys.exit()
else:
    pass

if cmd.template is not None:
    with open(cmd.template, 'w') as configtemplate:
        configtemplate.write(tupa_help.template_content)
    sys.exit("\n>>> Configuration template file (" + str(cmd.template) + ") written.\n>>> Exiting...\n")
else:
    pass


#####################################
# Parsing command line arguments
top_file      = cmd.top
traj_file     = cmd.traj
outdir        = cmd.outdir
configinput   = os.path.abspath(cmd.config)
if cmd.dumptime is None:
    dumptime = []
else:
    dumptime      = cmd.dumptime

n_cpus = int(cmd.ncpus)

#####################################
# Parsing configuration file and parameters
config = Configuration(configinput)


#####################################
# writing and printing run info
format.print_run_info(config, top_file, traj_file, outdir)

# creating file objects to write in them on the fly
out, outres, outangle, outprobe = format.open_outputs(outdir)
if config.mode == "bond":
    outproj, outalig, outrbond, outresbond = format.open_outputs_bond(outdir)


###############################################################################
###############################################################################
# From this point on, we will start to 1) process topology and trajectory;
# 2) setup some conditions to proceed with the calculation; 3) actually
# calculate the electric field in different modes and 4) write to the
# different output files.
###############################################################################
###############################################################################


#####################################
# Creating MDA Universe
print("\n########################################################")
if top_file != None and traj_file != None:
    try:
        u = mda.Universe(top_file, traj_file)
    except:
        raise ValueError(""">>> ERROR: Topology or Trajectory could not be loaded. Check your inputs!\n""")
else:
    raise ValueError(""">>> ERROR: Topology or Trajectory not found. Aren't you forgetting something?!\n""")


# Check whether trajectory file has box dimension information and redefine it
# (if it has been requested in the configuration file)
always_redefine_box_flag = check_box_vectors(u, config)
if always_redefine_box_flag:
    u.dimensions = config.boxdimensions

#####################################
# Sanity check1: sele_environment and probe_selection should NOT be empty
# Selecting the atoms that exert the electric field
try:
    env_raw_selection = u.select_atoms(config.sele_elecfield)
    if not env_raw_selection:
        # Sanity check1: environment selection can not be empty.
        raise ValueError(">>> ERROR: ENVIRONMENT selection (from sele_environment) is empty. Check your selection!\n")
except:
    # Sanity check1: environment selection must be parsed.
    raise ValueError(">>> ERROR: ENVIRONMENT selection (from sele_environment) cannot be parsed. Check your selection!\n")

probe_selection = create_probe_selection(u, config)


#####################################
# Sanity check2: atoms in probe_selection should NOT be in env_raw_selection
if config.mode == "atom" or config.mode == "bond":
    for atom in probe_selection.atoms:
        if atom in env_raw_selection.atoms:
            print(">>> WARNING: Probe atom (" + str(atom) + ") included in Environment selection! Moving on...!\n")


###############################################################################
# Looping over the trajectory
print("\n########################################################")
print("\n>>> Calculating Electric Field:")

results = Results() # Creating the Result Object carrying the analyzes results


for ts in u.trajectory:
    #####################################
    # Update environment selections and probe position
    refposition = update_probe_position(u, config, probe_selection)
    # We incorporate the solvent selection around the probe here for each mode
    enviroment_selection = update_environment(u, config, env_raw_selection, refposition)

	#####################################
    # Converting MDanalysis frames using defined dt in picoseconds
    time = config.dt * (int(ts.frame) + 1)
    print("Time (ps) = {0:<12d} -> (Probe position: {1})".format(time,refposition))

    # Write the probe coordinates
    lineprobe = format.fmt_probe_position(time, refposition)
    outprobe.write(lineprobe)

    # Dumping a specific frame if asked
    dump_coor_time(u, time, dumptime, enviroment_selection, outdir)

    #####################################
    # Using multiple processes to calculate atomic contribution to the
    # electric field. Store results as an array like:
    # [ [ atomid1   array([Xvec1  Yvec1  Zvec1]) ]
    #   [ atomid2   array([Xvec2  Yvec2  Zvec2]) ]
    #    ...
    #   [ atomidN   array([XvecN  YvecN  ZvecN]) ]   ]

    def calc_atom_contribution(atomid, refposition):

        atom = enviroment_selection.atoms[atomid]
        Ef_xyz = calc_ElectricField(atom, refposition)

        return [atomid,Ef_xyz]

    atomid_values = np.arange(len(enviroment_selection.atoms))
    run_per_atom = partial(calc_atom_contribution, refposition = refposition)

    # create a default thread pool
    with Pool(n_cpus) as worker_pool:
        pool_results = worker_pool.map(run_per_atom, atomid_values)

    per_atom_Ef = np.vstack(np.array(pool_results, dtype=object))[:,1] #ignore atomids
    totalEf  = np.sum(per_atom_Ef, axis=0)
    totalEfmag = mag(totalEf)

    # Keep the contributions for each frame so we can use later
    results.efield_timeseries[time] = np.append(totalEfmag, totalEf)

    # write Efield
    lineEfield  = format.fmt_efield_out_line(results.efield_timeseries[time], t=time)
    out.write(lineEfield)

    #####################################
    # Calculate information exclusive to BOND mode (Projection of Efield to bond)
    if config.mode == "bond":
        Efproj, proj_direction, angle_deg, rbond_vec = calculate_Efprojection(u, totalEf, config)
        Efprojmag = mag(Efproj)*proj_direction

        results.projection_timeseries[time] = np.append(Efprojmag, Efproj)

        # Write Efield projection
        lineproj = format.fmt_proj_out_line(results.projection_timeseries[time], t=time, optarr=[proj_direction])
        outproj.write(lineproj)

        # Calculate and write percentage of projection alignment
        aligned = alignment(Efprojmag,totalEfmag)*proj_direction

        results.projalignment[time] = [angle_deg, aligned]
        linealig = format.fmt_align_out_line([time, angle_deg, aligned])
        outalig.write(linealig)

        # Write rbond_info data
        bondline = format.fmt_rbond_info([time, rbond_vec, unit_vector(rbond_vec), mag(rbond_vec)])
        outrbond.write(bondline)


    #####################################
    # opening a temporary dictionary to hold the contribution of each residue
    # for each frame
    tmp_dict_res = {}
    for atomic_contribution in pool_results:
        atomid = atomic_contribution[0]
        contribution = atomic_contribution[1]

        atom = enviroment_selection.atoms[atomid]
        residue = str(atom.segid) + "_" + str(atom.resid)

        if atom in env_raw_selection.atoms:
            tmp_dict_res = add_to_dict(tmp_dict_res, residue, contribution)

    #####################################
    # Update dictionary with the contribution of each residue with the
    # values of Efield of each residue in this frame
    for r, contribution in tmp_dict_res.items():

        # account for residues containing 1 single atom
        if len(contribution.shape) > 1:
            resEf  = np.sum(contribution, axis=0)
        else:
            resEf = contribution
        resEfmag = mag(resEf)

        # start residue array
        resEf_array = np.append(resEfmag, resEf)
        # calculate the projection and alignment for each residue
        resEfproj = projection(resEf,totalEf) # resEf can have a higher magnitude than the total field.
        resEfprojmag = mag(resEfproj)
        resEfalignment = alignment(resEfprojmag,totalEfmag) # resEf can have a higher magnitude than the total field,
                                                              # which means % can be > 100%

        resEf_array = np.append(resEf_array,[resEfprojmag,resEfalignment])
        results.res_contribution_per_frame = add_to_dict(results.res_contribution_per_frame, r, resEf_array)


        # Calculate residue contribution to the projected efield
        if config.mode == "bond":
            resEfproj_bond = projection(resEf,rbond_vec)

            # Calculate projection direction (either 1 or -1)
            angle_deg = angle_between(resEfproj_bond,rbond_vec)
            proj_direction = np.cos(angle_deg)/abs(np.cos(angle_deg))

            resEfprojmag_bond = mag(resEfproj_bond)*proj_direction
            resEfalignment_bond = alignment(resEfprojmag_bond,totalEfmag)            

            if r not in results.resEFalignment_bond_per_frame.keys():
                new_res_array1 = np.append(time, [resEfprojmag_bond,resEfalignment_bond])
                results.resEFalignment_bond_per_frame[r]  = [new_res_array1]
            else:
                new_res_array1 = np.append(time, [resEfprojmag_bond,resEfalignment_bond])
                old_res_array  = results.resEFalignment_bond_per_frame[r]
                new_res_array2 = np.append(old_res_array, [new_res_array1], axis=0)
                results.resEFalignment_bond_per_frame[r]  = new_res_array2


###############################################################################
# Calculate the average contribution of each residue to the total field
for r, timeseries in results.res_contribution_per_frame.items():
    # timeseries[0] = resEfmag
    # timeseries[1,2,3] = resEfmag
    # timeseries[4] = resEfprojmag
    # timeseries[5] = resEfalignment
    if len(timeseries.shape) > 1: # to handle multiple frames (2D array)
        resEfmag       = timeseries[:,0].astype('float32')
        resEfprojmag   = timeseries[:,4].astype('float32')
        resEfalignment = timeseries[:,5].astype('float32')

        resEfmag_avg  = np.average(resEfmag)
        resEfmag_std  = np.std(resEfmag)
        resEfalig_avg = np.average(resEfalignment)
        resEfalig_std = np.std(resEfalignment)
    else:                         # to handle a single frame/pdb (1D array)
        resEfmag       = timeseries[0].astype('float32')
        resEfprojmag   = timeseries[4].astype('float32')
        resEfalignment = timeseries[5].astype('float32')

        resEfmag_avg  = np.average(resEfmag)
        resEfmag_std  = np.std(resEfmag)
        resEfalig_avg = np.average(resEfalignment)
        resEfalig_std = np.std(resEfalignment)

    # write contribution of each residue
    lineres = format.fmt_res_out_line([r, resEfmag_avg, resEfmag_std, resEfalig_avg, resEfalig_std])
    outres.write(lineres)

if config.mode == 'bond':
    for r, timeseries in results.resEFalignment_bond_per_frame.items():
        #r = r.partition('_')[2] # ignore resname and keep resid
        resEfprojmag_bond   = timeseries[:,1].astype('float32')
        resEfalignment_bond = timeseries[:,2].astype('float32')

        resEfmag_avg_bond  = np.average(resEfprojmag_bond)
        resEfmag_std_bond  = np.std(resEfprojmag_bond)
        resEfalig_avg_bond = np.average(resEfalignment_bond)
        resEfalig_std_bond = np.std(resEfalignment_bond)

        # write contribution of each residue
        lineresbond = format.fmt_res_out_line([r, resEfmag_avg_bond, resEfmag_std_bond, resEfalig_avg_bond, resEfalig_std_bond])
        outresbond.write(lineresbond)


###############################################################################
# Calculate average Efield vector and its angle to the field vector of each frame (Efield(t)).
# The ultimate goal here is to create a 3D standard deviation to be plotted with the field vector in pyTUPÃ.
avgfield        = np.average(list(results.efield_timeseries.values()), axis=0)
stdfield        = np.std(list(results.efield_timeseries.values()), axis=0)
avgfieldmag     = avgfield[0] # magnitude is the first array item


results.spatialdev = np.array([], dtype='float32')

for time,field in results.efield_timeseries.items():
    # Projection between Efield(t) and average Efield
    tmp_proj      = projection(field[1:],avgfield[1:]) #excluding the first item (magnitude)
    tmp_projmag   = mag(tmp_proj)
    # Alignment between Efield(t) and average Efield
    tmp_projalig     = alignment(tmp_projmag,avgfieldmag)
    #theta = Angle between Efield(t) and average Efield
    theta     = angle_between(field, avgfield)

    tmp_spacedev_array = np.array([theta, tmp_projmag, tmp_projalig])
    results.spatialdev = np.append(results.spatialdev, [tmp_spacedev_array])

    lineangle  = format.fmt_spatialdev_out_line(tmp_spacedev_array, t=time)
    outangle.write(lineangle)

# reshape the array to allow sum vertically
results.spatialdev = np.reshape(results.spatialdev, (-1,3))
# write average angle between Efield(t) and the average Efield.
avgangle,   avgprojmag, avgprojalig = np.average(results.spatialdev, axis=0)
stdevangle, stdprojmag, stdprojalig = np.std(results.spatialdev, axis=0)

outangle.write("#---#\n")
string_avgangle = format.fmt_spatialdev_out_line([avgangle,   avgprojmag, avgprojalig], lastline = True)
string_stdangle = format.fmt_spatialdev_out_line([stdevangle, stdprojmag, stdprojalig], lastline = True)
outangle.write("#AVG:       " + string_avgangle)
outangle.write("#STDEV:     " + string_stdangle)
outangle.close()

# write average Efield to the output file
avgmag, avgx,  avgy,  avgz   = avgfield
stdmag, stdx , stdy,  stdz   = stdfield

out.write("#---#\n")
string_avgmag = format.fmt_efield_out_line([avgmag, avgx, avgy, avgz], lastline = True)
string_stdmag = format.fmt_efield_out_line([stdmag, stdx, stdy, stdz], lastline = True)
out.write("#AVG:       " + string_avgmag)
out.write("#STDEV:     " + string_stdmag)

if config.mode == "bond":
    avgproj         = np.average(list(results.projection_timeseries.values()), axis=0)
    stdproj         = np.std(list(results.projection_timeseries.values()), axis=0)
    # write average ElecField_proj_onto_bond to the output file
    avgmag, avgx, avgy, avgz   = avgproj
    stdmag, stdx, stdy, stdz   = stdproj

    outproj.write("#---#\n")
    string_avgmag = format.fmt_proj_out_line([avgmag, avgx, avgy, avgz], lastline = True)
    string_stdmag = format.fmt_proj_out_line([stdmag, stdx, stdy, stdz], lastline = True)
    outproj.write("#AVG:       " + string_avgmag)
    outproj.write("#STDEV:     " + string_stdmag)

    # write average alignment and average angle between Efield(t) and Proj(t)
    avgangle,   avgalign   = np.average(list(results.projalignment.values()), axis=0)
    stdevangle, stdevalign = np.std(list(results.projalignment.values()), axis=0)

    outalig.write("#---#\n")
    string_avgalig = format.fmt_align_out_line([avgangle,   avgalign],   lastline = True)
    string_stdalig = format.fmt_align_out_line([stdevangle, stdevalign], lastline = True)
    outalig.write("#AVG:       " + string_avgalig)
    outalig.write("#STDEV:     " + string_stdalig)
    outalig.close()


###############################################################################
# Close output files
outres.close()
out.close()
if config.mode == "bond":
    outproj.close()
    outrbond.close()

stop = timeit.default_timer()
print("\n>>> Calculation complete!")
print('>>> Runtime: {} sec'.format(round(stop - start,2)))
