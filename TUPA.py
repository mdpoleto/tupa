#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#Marcelo D. Poleto
#FEB 2022

import sys, os, argparse, timeit
import MDAnalysis as mda
import numpy as np
import configparser as cp
import warnings
sys.dont_write_bytecode = True
warnings.filterwarnings("ignore", message="Found no information for attr:")
warnings.filterwarnings("ignore", message="Found missing chainIDs")
np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

import tupa_help

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

cmd = ap.parse_args()
start = timeit.default_timer()
#print(tupa_help.header)

if cmd.help is True:
	ap.print_help()
	print(tupa_help.help)
	sys.exit()
else:
	pass

if cmd.template is not None:
	with open(cmd.template, 'w') as configtemplate:
		configtemplate.write(tupa_help.template_content)
	print("\n>>> Configuration template file (" + str(cmd.template) + ") written.\n>>> Exiting...\n")
	sys.exit()
else:
	pass


top_file      = cmd.top
traj_file     = cmd.traj
outdir        = cmd.outdir
configinput   = cmd.config
if cmd.dumptime is None:
	dumptime = []
else:
	dumptime      = cmd.dumptime
###############################################################################
# Stablishing default config values
config = cp.ConfigParser(allow_no_value=True, inline_comment_prefixes="#")

###############################################################################
# Reading configuration file and assuming absences
if configinput == None:
	print("\n>>> Configuration file is a necessary input. Use -template to generate one.\n>>> Exiting...\n")
	sys.exit()
else:
	config.read(configinput)

try:
	sele_elecfield = str(config['Environment Selection']["sele_environment"].strip('"'))
except Exception as e1:
	sys.exit("\n>>> ERROR: sele_environment must be defined with a valid MDanalysis selection!\n")

try:
	mode = str(config['Probe Selection']["mode"].strip('"')).lower()
	if mode == "atom":
		try:
			selatom = str(config['Probe Selection']["selatom"].strip('"'))
		except Exception as e2:
			sys.exit("""\n>>> ERROR: in "ATOM" mode, selatom must be defined!\n""")
	elif mode == "bond":
		try:
			selbond1 = str(config['Probe Selection']["selbond1"].strip('"'))
			selbond2 = str(config['Probe Selection']["selbond2"].strip('"'))
		except Exception as e2:
			sys.exit("""\n>>> ERROR: in "BOND" mode, both selbond1 and selbond2 must be defined!\n""")
	elif mode == "coordinate":
		try:
			tmp_probecoordinate = str(config['Probe Selection']["probecoordinate"]).strip('[]').split(",")
			probecoordinate = [float(item) for item in tmp_probecoordinate]
		except Exception as e2:
			sys.exit("""\n>>> ERROR: in "COORDINATE" mode, a list of coordinates [X,Y,Z] must be provided!\n""")
	elif mode == "list":
		try:
			file_of_coordinates = str(config['Probe Selection']["file_of_coordinates"].strip('"'))
		except Exception as e2:
			sys.exit("""\n>>> ERROR: in "LIST" mode, "file_of_coordinates" must be defined!\n""")
	else:
		sys.exit("""\n>>> ERROR: "mode" must be defined as "ATOM", "BOND", "COORDINATE" or "LIST"!\n""")

except Exception as e2:
	sys.exit("""\n>>> ERROR: "mode" must be defined as "ATOM", "BOND", "COORDINATE" or "LIST"!\n""")


if mode == "coordinate" or mode == "list":
	try:
		remove_self = str(config['Probe Selection']["remove_self"]).lower()
		if remove_self  == "true":
			remove_self = True
			try:
				remove_cutoff = float(config['Probe Selection']["remove_cutoff"])
			except Exception as e3:
				print("""\n>>> WARNING: "remove_cutoff" was not provided! Using the default value (1 Angstrom)!...\n""")
				remove_cutoff = 1
		else:
			remove_self = False
			remove_cutoff = 0
	except:
		remove_self = False
		remove_cutoff = 0
else:
	remove_self = False
	remove_cutoff = 0


try:
	include_solvent = str(config['Solvent']["include_solvent"].strip('"')).lower()
	if include_solvent == "true":
		include_solvent   = True
		try:
			solvent_selection = str(config['Solvent']["solvent_selection"].strip('"'))
		except Exception as e4:
			sys.exit("""\n>>> ERROR: solvent_selection must be a valid MDanalysis selection!\n""")
		try:
			solvent_cutoff    = float(config['Solvent']["solvent_cutoff"])
		except Exception as e5:
			sys.exit("""\n>>> ERROR: solvent_cutoff must be a number!\n""")
	else:
		include_solvent   = False
except:
	include_solvent = False


try:
	dt = int(config['Time']["dt"])
except:
	dt = 1

try:
	redefine_box = str(config['Box Info']["redefine_box"].strip('"')).lower()
	if redefine_box == "true":
		redefine_box = True
		try:
			tmp_boxdim = str(config['Box Info']["boxdimensions"].strip('[]')).split(",")
			boxdimensions = [float(item) for item in tmp_boxdim]
			if len(boxdimensions) != 6:
				sys.exit("""\n>>> ERROR: boxdimensions should contain [a,b,c, alpha, beta, gamma]!\n""")
		except Exception as e6:
			sys.exit("""\n>>> ERROR: To redifine the box, boxdimensions must be provided!\n""")
	else:
		redefine_box = False
		boxdimensions = None
except:
	redefine_box = False
	boxdimensions = None

###############################################################################
# Being verbose about parameters chosen
print("########################################################")
print(">>> Parameters used to run TUPÃ:")
print()
print('Topology file       = ' + top_file)
print('Trajectory file     = ' + traj_file)
print()
print('[Environment Selection]')
print('sele_environment    = {}'.format(sele_elecfield))
print()
print("[Probe Selection]")
print('mode                = {}'.format(mode))
probestr = ""

if mode == "atom":
	str1 = 'selatom             = {}'.format(selatom)
	print(str1)
	probestr += str1
elif mode == "bond":
	str1 = 'selbond1            = {}'.format(selbond1)
	str2 = 'selbond2            = {}'.format(selbond2)
	print(str1)
	print(str2)
	probestr += str1 +'\n'+str2
elif mode == "coordinate":
	str1 = 'probecoordinate     = {}'.format(probecoordinate)
	print(str1)
	probestr += str1
elif mode == "list":
	str1 = 'file_of_coordinates = {}'.format(file_of_coordinates)
	print(str1)
	probestr += str1+'\n'
	if remove_self == True:
		str2 = 'remove_self         = {}'.format(remove_self)
		str3 = 'remove_cutoff       = {}'.format(remove_cutoff)
		print(str2)
		print(str3)
		probestr += str2 +'\n' + str3
print()
print("[Solvent]")
print('include_solvent     = {}'.format(include_solvent))
solvstr = ""
if include_solvent == True:
	str1 = 'solvent_selection   = {}'.format(solvent_selection)
	str2 = 'solvent_cutoff      = {}'.format(solvent_cutoff)
	print(str1)
	print(str2)
	solvstr += str1 +"\n"+str2
print()
print("[Time]")
print('dt                  = {}'.format(dt))
print()
print("[Box Info]")
print('redefine_box        = {}'.format(redefine_box))
boxstr = ""
if redefine_box == True:
	str1 = 'boxdimensions       = {}'.format(boxdimensions)
	print(str1)
	boxstr += str1

###############################################################################
def mag(vector):
	""" Returns the magnitude of the vector.  """
	mag = np.linalg.norm(vector)
	return mag

def unit_vector(vector):
	""" Returns the unit vector of the vector.  """
	return vector / np.linalg.norm(vector)

def alignment(value1,value2):
	# V2 is the total reference
	aligned = abs(value1/value2)*100
	return aligned

def projection(v1,v2):
	#projection_u_on_v = (np.dot(u, v)/np.dot(v, v))*v
	proj = v2*(np.dot(v1, v2)/np.dot(v2,v2))
	return proj

def angle_between(v1, v2):

	v1_u = unit_vector(v1)
	v2_u = unit_vector(v2)

	angle_rad = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
	return angle_rad

def calc_EletricField(atom,refposition):

	# Set Epsilon to your medium (in C**2/N*(m**2))
	Epsilon = 8.8541878128*(10**(-12))
	k = 1/(4*np.pi*Epsilon)  # 8987551792.261173  (N*m**2)/C**2

	# Electric field = N/C = V/m

	refX = refposition[0]*(10**(-10)) #convert from Angstrom to meter
	refY = refposition[1]*(10**(-10)) #convert from Angstrom to meter
	refZ = refposition[2]*(10**(-10)) #convert from Angstrom to meter

	atomX, atomY, atomZ = atom.position
	atomX = atomX*(10**(-10)) #convert from Angstrom to meter
	atomY = atomY*(10**(-10)) #convert from Angstrom to meter
	atomZ = atomZ*(10**(-10)) #convert from Angstrom to meter

	r2 = [refX, refY, refZ]
	r1 = [atomX, atomY, atomZ]

	rvec = np.array(r2) - np.array(r1)
	rmag = mag(rvec)
	rhat = unit_vector(rvec)

	# convert charge (elementary charge unit to Coulomb)
	charge = round(atom.charge, 6)   # round to 6 digits so the net charge is closer to 0
	# 1 elementary charge unit = 1.60217733*(10**(-19))
	charge = charge*(1.60217733*(10**(-19))) # convert to Coulomb

	# Calculate Electric Field
	# E = k*Q/r2; charge = coulomb ; r = meter, then: Ef = N/C
	Ef = rhat*(k*charge/rmag**2)
	Ef = Ef*(10**(-8)) #Convert to MegaVolt per centimeter (usual unit in literature)
	# 1 MV/cm equals 10**8 V/m (the SI unit), 1.944*10**(−4) atomic units, or 0.0480 kcal/(mol*Debye)
	# 10 MV/cm = 1 V/nm = 0.1 V/Angstrom

	return Ef

def pack_around(atom_group, center, boxinfo):

	tmp_group = atom_group # we will update tmp_group below
	atom_group.pack_into_box()
	positions = atom_group.positions.copy() # working with a tmp copy of positions
	sub = positions - center

	# Get the box for the current frame
	box = boxinfo[:3]
	culprits = np.where(np.sqrt(sub**2) > box / 2)
	# Actually translate the coordinates.
	positions[culprits] -= (u.dimensions[culprits[1]]
	                        * np.sign(sub[culprits]))
	tmp_group.positions = positions

	return tmp_group
###############################################################################
# Opening output files

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
""".format(top = top_file,
           traj = traj_file,
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


out = open(os.path.join(outdir,"ElecField.dat"),'w')
out.write("""@    title "Electric Field"\n@    xaxis  label "Time (ps)"\n@    yaxis  label "MV/cm"\n""")
out.write("#time      Magnitude          Efield_X            Efield_Y            Efield_Z\n")
out.write("@type xy\n")

outres = open(os.path.join(outdir,"ElecField_per_residue.dat"),'w')
outres.write("""@    title "Magnitude of Electric Field"\n@    xaxis  label "Residues"\n@    yaxis  label "MV/cm"\n""")
outres.write("#time     Efield_per_residue  Std.Deviation       Alignment_percent   Alignment_Stdev\n")
outres.write("@type xydy\n")

outangle = open(os.path.join(outdir,"Spatial_Deviation.dat"), "w")
outangle.write("#time   Angle(Efield(t), avg_field)   Projection(Efield(t), avg_field)   Alignment(Efield(t), avg_field)\n")

outprobe = open(os.path.join(outdir,"probe_info.dat"), "w")
outprobe.write("#time     probe_coordinates\n")

if mode == "bond":
	outproj = open(os.path.join(outdir,"ElecField_proj_onto_bond.dat"),'w')
	outproj.write("""@    title "Efield_Projection"\n@    xaxis  label "Time (ps)"\n@    yaxis  label "MV/cm"\n""")
	outproj.write("#time      Magnitude          Efield_X            Efield_Y            Efield_Z          Direction\n")
	outproj.write("@type xy\n")

	outalig = open(os.path.join(outdir,"ElecField_alignment.dat"),'w')
	outalig.write("""@    title "Efield_Projection"\n@    xaxis  label "Time (ps)"\n@    yaxis  label "Alignment_percentage"\n""")
	outalig.write("#time     Alignment_percent   Angle\n")
	outalig.write("@type xy\n")

	outrbond = open(os.path.join(outdir,"bond_axis_info.dat"), "w")
	outrbond.write("#time     rbondvec_in_meters                                          rbond_magnitude               rbond_unit_vector\n")

###############################################################################
print("\n########################################################")
# Creating Universe and making selections
if top_file != None and traj_file != None:
	try:
		u = mda.Universe(top_file, traj_file)
	except Exception as e:
		sys.exit("""\n>>> ERROR: Topology or Trajectory could not be loaded. Check your inputs!\n""")

else:
	sys.exit("""\n>>> ERROR: Topology or Trajectory not found. Aren't you forgetting something?!\n""")

###############################################################################
print("\n########################################################")
# Check whether trajectory file has box dimension information and redefine if requested
if u.dimensions is not None:
	boxangles = u.dimensions[3:]
	checkangles = [i for i in boxangles if float(i) != 90.]
	if len(checkangles) == 0:
		if u.dimensions[0] == 1 and u.dimensions[1] == 1 and u.dimensions[2] == 1:
			if redefine_box == True:
				print("\n>>> Redefining box dimensions to:", boxdimensions)
				always_redefine_box_flag = True
			else:
				sys.exit("""\n>>> ERROR: Your trajectory does not contain information regarding box size. Provide them in the configuration file!\n""")
		elif u.dimensions[0] == 0 and u.dimensions[1] == 0 and u.dimensions[2] == 0:
				print("\n>>> Redefining box dimensions to: ", boxdimensions)
				always_redefine_box_flag = True
		else:
			if redefine_box == True:
				print("\n>>> Redefining box dimensions to:", boxdimensions)
				always_redefine_box_flag = True
			else:
				always_redefine_box_flag = False
	else:
		sys.exit("""\n>>> WARNING: Your box does not seem to be orthorhombic. We are currently working to support non-rectangular boxes.\n""")
		always_redefine_box_flag = False
elif u.dimensions is None and redefine_box is True:
	always_redefine_box_flag = True
else:
	sys.exit("""\n>>> ERROR: Your trajectory does not contain information regarding box size. Provide them in the configuration file!\n""")

###############################################################################
# This is us being very verbose so people actually know what is happening
if mode == "atom":
	try:
		probe_selection = u.select_atoms(selatom, updating=True)
	except Exception as e:
		###############################################################################
		# Sanity check1: environment selection can not be empty.
		sys.exit("\n>>> ERROR: PROBE selection (from selatom) is empty. Check your selection!\n")

	if len(probe_selection) > 1:
		print("\n>>> Running in ATOM mode!")
		print(">>> Probe atoms (using center of geometry) = "+ str(probe_selection.atoms) + "\n")
	else:
		print("\n>>> Running in ATOM mode!")
		print(">>> Probe atom = "+ str(probe_selection.atoms[0].name) +"-"+ str(probe_selection.resnames[0]) + str(probe_selection.resnums[0]) + "\n")
elif mode == "bond":
	try:
		bond1 = u.select_atoms(selbond1, updating=True)
	except Exception as e:
		###############################################################################
		# Sanity check1: PROBE selection can not be empty.
		sys.exit("\n>>> ERROR: PROBE selection (from selbond1) is empty. Check your selection!\n")
	try:
		bond2 = u.select_atoms(selbond2, updating=True)
	except Exception as e:
		###############################################################################
		# Sanity check1: PROBE selection can not be empty.
		sys.exit("\n>>> ERROR: PROBE selection (from selbond2) is empty. Check your selection!\n")

	probe_selection = bond1 + bond2
	print("\n>>> Running in BOND mode!")
	print(">>> Probe atoms = " + str(probe_selection.atoms))
	print(">>> Bond axis direction = " + str(bond1.atoms[0].name) +"-"+ str(bond1.resnames[0]) + str(bond1.resnums[0]) + " --> " + str(bond2.atoms[0].name) +"-"+ str(bond2.resnames[0]) + str(bond2.resnums[0]) + "\n")
elif mode == "coordinate":
	probe_selection = probecoordinate
	print("\n>>> Running in COORDINATE mode!")
	print(">>> Probe XYZ position = " +  str(probe_selection))
	if remove_self == True:
		print(">>> Removing self contribution within a radial cutoff of " +  str(remove_cutoff) + " Angstroms\n")
elif mode == "list":
	print("\n>>> Running in LIST mode!")
	# If mode == LIST, parse the file of coordinates
	loc = open(file_of_coordinates, 'r')
	probe_selection = []
	while True:
		coorline = loc.readline()
		if not coorline: break
		if coorline[0] != "#" and coorline[0] != "@":
			if len(coorline.split()) == 3:
				try:
					listcoorX = float(coorline.split()[0])
					listcoorY = float(coorline.split()[1])
					listcoorZ = float(coorline.split()[2])
				except:
					sys.exit("""\n>>> ERROR: File of coordinates could not be parsed! Expecting floats or integers (X Y Z). Check your inputs!\n""")
				tmpcoorlist = [listcoorX, listcoorY, listcoorZ]
				probe_selection.append(tmpcoorlist)
			else:
				sys.exit("""\n>>> ERROR: File of coordinates could not be parsed! Expecting 3 columns (X Y Z). Check your inputs!\n""")

	# Sanity check: list of coordinates should have the same length as the trajectory
	if len(probe_selection) != len(u.trajectory):
		sys.exit("""\n>>> ERROR: File of coordinates and trajectory have different lengths (""" + str(len(probe_selection)) + """ lines and """ + str(len(u.trajectory)) + """ frames)!\n""")

	print(">>> Probe XYZ positions = " + str(len(probe_selection)) + """ XYZ coordinates""")
	if remove_self == True:
		print(">>> Removing self contribution within a radial cutoff of " +  str(remove_cutoff) + " Angstroms\n")

else:
	sys.exit("\n>>> MODE should be 'bond', 'atom', 'coordinate' or 'list'!")


###############################################################################
# Selecting the atoms that exert the electric field
try:
	elecfield_selection = u.select_atoms(sele_elecfield)
except Exception as e:
	###############################################################################
	# Sanity check1: environment selection can not be empty.
	sys.exit("\n>>> ERROR: ENVIRONMENT selection (from sele_environment) is empty. Check your selection!\n")

###############################################################################
# Sanity check2: atoms in probe_selection should NOT be in elecfield_selection
# if tmprefposition == probe_selection.center_of_geometry(). Checking that...
if mode == "atom" or mode == "bond":
	for atom in probe_selection.atoms:
		if atom in elecfield_selection.atoms:
			print(">>> WARNING: Probe atom (" + str(atom) + ") included in Environment selection! Moving on...!\n")

###############################################################################
# Verbose output for solvent inclusion in calculation. Selection is done within the trajectory loop
if include_solvent == True:
	print(">>> Including solvent molecules within a radius of " + str(solvent_cutoff) + " Angstroms in the calculation...")
else:
	print(">>> Excluding solvent molecules from the calculation...")

###############################################################################
# Opening a dictionary to hold residues and their Efield contribution
dict_res_total = {}

for atom in elecfield_selection.atoms:
	residue = str(atom.resname) + "_" + str(atom.resid)
	if residue not in dict_res_total.keys():
		dict_res_total[residue] = [ [],[],[],[],[] ]

###############################################################################
print("\n########################################################")
print("\n>>> Calculating Electric Field:")

efield_total     = {}
projection_total = {}
align_list       = []
angle_list       = []
efieldmag_list   = []
projmag_list     = []


for ts in u.trajectory[0: len(u.trajectory):]:
	# if we detected that box needs dimensions, we redefine it on each frame here
	if always_redefine_box_flag == True:
		ts.dimensions = boxdimensions
	else:
		pass

	########################################################
	# Update environment selections and probe position

	if len(probe_selection) == 0:
		sys.exit("\n>>> ERROR: PROBE selection is empty! Check your selection!\n")


	if mode == "atom":
		if len(probe_selection) > 1:
			refposition = probe_selection.center_of_geometry()
		else:
			refposition = probe_selection.atoms[0].position

	elif mode == "bond":
		position1 = bond1.atoms[0].position
		position2 = bond2.atoms[0].position

		rbond_vec = (position2 - position1) # axis from ref1 to ref2
		rbond_vec = rbond_vec*(10**(-10)) # convert from Angstrom to meter
		rbond_hat = unit_vector(rbond_vec)

		refposition = (position1 + position2)/2 # midway for both atoms

	elif mode == "coordinate":
		refposition = probe_selection

	elif mode == "list":
		refposition = probe_selection[ts.frame] # update probe position from the list

	# We incorporate the solvent selection around the probe here for each mode
	if include_solvent == True:
		if mode == "atom":
			tmp_selection = u.select_atoms("(around " + str(solvent_cutoff) + " " + selatom + ") and (" + solvent_selection + ")", periodic=True)
		elif mode == "bond":
			tmp_selection = u.select_atoms("(around " + str(solvent_cutoff) + " (" + selbond1 + " or " + selbond2 +")) and (" + solvent_selection + ")", periodic=True)
		elif mode == "coordinate":
			tmp_selection = u.select_atoms("(point " + str(refposition[0]) + " " + str(refposition[1]) + " " + str(refposition[2]) + " " + str(solvent_cutoff) + ") and (" + solvent_selection + ")", periodic=True)
		elif mode == "list":
			tmp_selection = u.select_atoms("(point " + str(refposition[0]) + " " + str(refposition[1]) + " " + str(refposition[2]) + " " + str(solvent_cutoff) + ") and (" + solvent_selection + ")", periodic=True)

		# translate molecules within selection that are beyond the PBC and incorporate into the environment
		tmp_selection = pack_around(tmp_selection, refposition, ts.dimensions)
		enviroment_selection = elecfield_selection + tmp_selection
	else:
		enviroment_selection = elecfield_selection

	########################################################

	# Converting MDanalysis frames using defined dt
	frame = int(ts.frame) + 1
	time = "{:.0f}".format(frame*dt)
	print("Time = " + str(time) + " ps... (Probe position: " + str(refposition) + ")")

	# Write the probe coordinates
	listprobe = "[" + str(refposition[0]) + "," + str(refposition[1]) + "," + str(refposition[2]) + "]"
	lineprobe = str(time).ljust(10,' ') +  listprobe.ljust(60,' ') + "\n"
	outprobe.write(lineprobe)

	# Dumping a specific frame if asked
	dumptime = np.array(dumptime, dtype=float)
	if float(time) in dumptime:
		framefile = os.path.join(outdir,"frame_" + time + "ps.pdb")
		envfile   = os.path.join(outdir,"environment_" + time + "ps.pdb")
		print("   >>> Dumping frame (Time = " + str(time) + " ps)! Check " + framefile)
		dump_sel = u.select_atoms("all")
		dump_sel.write(framefile)
		dump_sel2 = enviroment_selection.select_atoms("all")
		dump_sel2.write(envfile)

	# Evaluate self_contribution removal
	# selects all atoms from environment within a cutoff of a point in space (point X Y Z cutoff)
	# Such approach is only available to COORDINATE mode. ATOM and BOND modes can achieve the same behavior
	# by adjusting the environment selection accordingly.
	self_contribution = enviroment_selection.select_atoms("point  " + str(refposition[0]) + " " + str(refposition[1]) + " " + str(refposition[2]) + "  " + str(remove_cutoff) + "  ", periodic=True)

	if len(self_contribution.atoms) > 0:
		if remove_self == True:
			print(""">>> WARNING! Removing self contribution of: """ + str(self_contribution.atoms))
			# Remove self_contribution for this frame
			enviroment_selection = enviroment_selection - self_contribution
		else:
			print(""">>> WARNING! Some atoms are closed than """ + str(remove_cutoff) + """ A : """ + str(self_contribution.atoms))

	########################################################
	# opening a temporary dictionary to hold the contribution of each residue for each frame
	dict_res_tmp = {}

	xfield_list_frame = []
	yfield_list_frame = []
	zfield_list_frame = []

	# Iterating over all atoms in environment selection to get their contributions
	for atom in enviroment_selection.atoms:
		residue = str(atom.resname) + "_" + str(atom.resid)
		#print(residue)
		if residue not in dict_res_tmp.keys() and residue in dict_res_total.keys():
			dict_res_tmp[residue] = [[],[],[]]
		else:
			pass

		Ef_xyz = calc_EletricField(atom, refposition)

		Efx, Efy, Efz = Ef_xyz
		xfield_list_frame.append(Efx)
		yfield_list_frame.append(Efy)
		zfield_list_frame.append(Efz)

		# upload keys (residues) in dict_res_tmp to include the contribution of
		# each atom in a residue. We will update the dictionary for each residue
		# in this frame on line 487
		if residue in dict_res_total.keys():
			res_tmp_x, res_tmp_y, res_tmp_z = dict_res_tmp[residue]
			res_tmp_x.append(Efx)
			res_tmp_y.append(Efy)
			res_tmp_z.append(Efz)
			dict_res_tmp[residue] = [res_tmp_x, res_tmp_y, res_tmp_z]
		else:
			pass

	# Sum the contributions of all atoms in each axis to get the resultant Efield
	totalEfx = np.sum(xfield_list_frame)
	totalEfy = np.sum(yfield_list_frame)
	totalEfz = np.sum(zfield_list_frame)
	totalEf  = np.array([totalEfx, totalEfy, totalEfz])
	totalEfmag = mag(totalEf)

	# Keep the contributions for each frame so we can use later
	efield_total[time] = [totalEfx, totalEfy, totalEfz]
	# Keep the magnitudes so we can use later
	efieldmag_list.append(totalEfmag)

	########################################################
	# Write information exclusive to BOND mode
	if mode == "bond":
		# Calculate efield projection
		Efprojection = projection(totalEf,rbond_vec)
		Efproj_x, Efproj_y, Efproj_z = Efprojection
		projection_total[time] = [Efproj_x, Efproj_y, Efproj_z]

		# Calculate projection direction (either 1 or -1)
		angle_rad = angle_between(totalEf,rbond_vec)
		proj_direction = np.cos(angle_rad)/abs(np.cos(angle_rad))
		# Update the Efield projection sign depending on its direction (multiplying by 1 or -1)
		#totalEfmag = totalEfmag*proj_direction
		Efprojectionmag = mag(Efprojection)*proj_direction
		# Keep the magnitudes so we can use later
		projmag_list.append(Efprojectionmag)

		# write Efield
		lineEfield  = str(time).ljust(10,' ') + str("{:.6f}".format(totalEfmag)).ljust(20,' ') + str("{:.6f}".format(totalEfx)).ljust(20,' ') + str("{:.6f}".format(totalEfy)).ljust(20,' ') + str("{:.6f}".format(totalEfz)).ljust(20,' ') + "\n"
		out.write(lineEfield)
		# Write Efield projection
		lineproj = str(time).ljust(10,' ') + str("{:.6f}".format(Efprojectionmag)).ljust(20,' ') + str("{:.6f}".format(Efproj_x)).ljust(20,' ') + str("{:.6f}".format(Efproj_y)).ljust(20,' ') + str("{:.6f}".format(Efproj_z)).ljust(20,' ') + str(proj_direction) + "\n"
		outproj.write(lineproj)

		# Calculate and write percentage of projection alignment
		aligned = alignment(Efprojectionmag,totalEfmag)*proj_direction
		align_list.append(aligned)
		angle_list.append(np.degrees(angle_rad))
		linealig = str(time).ljust(10,' ') + str("{:.6f}".format(aligned)).ljust(20,' ') + str("{:.6f}".format(np.degrees(angle_rad))).ljust(20,' ') + "\n"
		outalig.write(linealig)

		# Write rbond_info data
		rvx, rvy, rvz = rbond_vec
		rhx, rhy, rhz = rbond_hat
		list1 = "[" + str(rvx) + "," + str(rvy) + "," + str(rvz) + "]"
		list2 = "[" + str(rhx) + "," + str(rhy) + "," + str(rhz) + "]"
		bondline = str(time).ljust(10,' ') +  list1.ljust(60,' ') +  str(mag(rbond_vec)).ljust(20,' ') + list2 + "\n"
		outrbond.write(bondline)

	else:
		# write Efield
		lineEfield  = str(time).ljust(10,' ') + str("{:.6f}".format(totalEfmag)).ljust(20,' ') + str("{:.6f}".format(totalEfx)).ljust(20,' ') + str("{:.6f}".format(totalEfy)).ljust(20,' ') + str("{:.6f}".format(totalEfz)).ljust(20,' ') + "\n"
		out.write(lineEfield)

	########################################################
	# Update dictionary with the contribution of each residue with the
	# values of Efield of each residue in this frame
	for r, comp in dict_res_tmp.items():
		res_efield_x, res_efield_y, res_efield_z = comp
		# values for the entire residue in THIS FRAME
		resEfx = np.sum(res_efield_x)
		resEfy = np.sum(res_efield_y)
		resEfz = np.sum(res_efield_z)
		resEf  = np.asarray([resEfx, resEfy, resEfz], dtype=object)
		resEfmag = mag(resEf)
		# calculate the projection and alignment for each residue
		if mode == "bond":
			resEfproj = projection(resEf,rbond_vec)  # resEf can have a higher magnitude than the total field.
		else:
			resEfproj = projection(resEf,totalEf)
		resEfprojmag = mag(resEfproj)
		resEfalignment = alignment(resEfprojmag,totalEfmag) # resEf can have a higher magnitude than the total field,
		                                                      # which means % can be > 100%

		# Update the total dictionary with values of THIS FRAME
		resEfx_tmp, resEfy_tmp, resEfz_tmp, resEfmag_tmp, res_alignment_tmp = dict_res_total[r]
		resEfx_tmp.append(resEfx)
		resEfy_tmp.append(resEfy)
		resEfz_tmp.append(resEfz)
		resEfmag_tmp.append(resEfmag)
		res_alignment_tmp.append(resEfalignment)
		dict_res_total[r] = [resEfx_tmp, resEfy_tmp, resEfz_tmp, resEfmag_tmp, res_alignment_tmp]

###############################################################################
# Calculate average Efield vector and its angle to the field vector of each frame (Efield(t)).
# The ultimate goal here is to create a 3D standard deviation to be plotted with the field vector in pyTUPÃ.
avgfield        = np.average(list(efield_total.values()), axis=0)
stdfield        = np.std(list(efield_total.values()), axis=0)
avgfieldmag     = mag(avgfield)

avgangle_list   = []
avgprojmag_list    = []
avgprojalig_list   = []


for time,field in efield_total.items():
	#Angle between Efield(t) and average Efield
	angle     = np.degrees(angle_between(field, avgfield))
	avgangle_list.append(angle)
	# Projection between Efield(t) and average Efield
	tmp_proj      = projection(field,avgfield)
	tmp_projmag   = mag(tmp_proj)
	avgprojmag_list.append(tmp_projmag)
	# Alignment between Efield(t) and average Efield
	tmp_projalig     = alignment(tmp_projmag,avgfieldmag)
	avgprojalig_list.append(tmp_projalig)

	lineangle  = str(time).ljust(10,' ') + str("{:.6f}".format(angle)).ljust(20,' ') + str("{:.6f}".format(tmp_projmag)).ljust(20,' ') + str("{:.6f}".format(tmp_projalig)).ljust(20,' ') + "\n"
	outangle.write(lineangle)

# write average angle between Efield(t) and the average Efield.
avgangle    = np.average(avgangle_list)
stdevangle  = np.std(avgangle_list)
avgprojmag  = np.average(avgprojmag_list)
stdprojmag  = np.std(avgprojmag_list)
avgprojalig = np.average(avgprojalig_list)
avgprojalig = np.std(avgprojalig_list)

outangle.write("#---#\n")
outangle.write("#AVG:   " + str("{:.2f}".format(avgangle)).rjust(6,' ') + " +- " + str("{:.2f}".format(stdevangle)).ljust(6,' ') + "      " + str("{:.2f}".format(avgprojmag)).rjust(6,' ') + " +- " + str("{:.2f}".format(stdprojmag)).ljust(6,' ') + "   " + str("{:.2f}".format(avgprojalig)).rjust(6,' ') + " +- " + str("{:.2f}".format(avgprojalig)).ljust(6,' '))
outangle.close()

# write average Efield to the output file
avgx,  avgy,  avgz   = avgfield
stdx , stdy,  stdz   = stdfield
avgmag = np.average(efieldmag_list)
stdmag = np.std(efieldmag_list)

out.write("#---#\n")
out.write("#AVG:     " + str("{:.6f}".format(avgmag)).ljust(20,' ') + str("{:.6f}".format(avgx)).ljust(20,' ') + str("{:.6f}".format(avgy)).ljust(20,' ') + str("{:.6f}".format(avgz)).ljust(20,' ') + "\n")
out.write("#STDEV:   " + str("{:.6f}".format(stdmag)).ljust(20,' ') + str("{:.6f}".format(stdx)).ljust(20,' ') + str("{:.6f}".format(stdy)).ljust(20,' ') + str("{:.6f}".format(stdz)).ljust(20,' ') + "\n")

if mode == "bond":
	avgproj         = np.average(list(projection_total.values()), axis=0)
	stdproj         = np.std(list(projection_total.values()), axis=0)
	# write average ElecField_proj_onto_bond to the output file
	avgx, avgy, avgz   = avgproj
	stdx, stdy, stdz   = stdproj
	avgmag = np.average(projmag_list)
	stdmag = np.std(projmag_list)

	outproj.write("#---#\n")
	outproj.write("#AVG:     " + str("{:.6f}".format(avgmag)).ljust(20,' ') + str("{:.6f}".format(avgx)).ljust(20,' ') + str("{:.6f}".format(avgy)).ljust(20,' ') + str("{:.6f}".format(avgz)).ljust(20,' ') + "\n")
	outproj.write("#STDEV:   " + str("{:.6f}".format(stdmag)).ljust(20,' ') + str("{:.6f}".format(stdx)).ljust(20,' ') + str("{:.6f}".format(stdy)).ljust(20,' ') + str("{:.6f}".format(stdz)).ljust(20,' ') + "\n")


	# write average alignment and average angle between Efield(t) and Proj(t)
	avgalign   = np.average(align_list)
	stdevalign = np.std(align_list)
	avgangle   = np.average(angle_list)
	stdevangle = np.std(angle_list)

	outalig.write("#---#\n")
	outalig.write("#AVG:     " + str("{:.2f}".format(avgalign)).rjust(6,' ') + " +- " + str("{:.2f}".format(stdevalign)).ljust(6,' ') + "    " + str("{:.2f}".format(avgangle)).rjust(6,' ') + " +- " + str("{:.2f}".format(stdevangle)).ljust(6,' '))
	outalig.close()
###############################################################################
# Calculate the average contribution of each residue to the total field
for r, comp in dict_res_total.items():
	r = r.partition('_')[2] # ignore resname and keep resid
	res_efield_x, res_efield_y, res_efield_z, resEfmag, resEfaligned = comp
	resEfmag_avg  = np.average(resEfmag)
	resEfmag_std  = np.std(resEfmag)
	resEfalig_avg = np.average(resEfaligned)
	resEfalig_std = np.std(resEfaligned)

	# write contribution of each residue
	lineres = str(r).ljust(10,' ') + str("{:.6f}".format(resEfmag_avg)).ljust(20,' ') + str("{:.6f}".format(resEfmag_std)).ljust(20,' ') + str("{:.6f}".format(resEfalig_avg)).ljust(20,' ') + str("{:.6f}".format(resEfalig_std)).ljust(20,' ') + "\n"
	outres.write(lineres)

###############################################################################
# Close output files
outres.close()
out.close()
if mode == "bond":
	outproj.close()
	outrbond.close()

stop = timeit.default_timer()
print("\n>>> Calculation complete!")
print('>>> Runtime: {} sec'.format(round(stop - start,2)))
