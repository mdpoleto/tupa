#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#Marcelo D. Poleto
#APRIL 2022

import sys, os, argparse, timeit
import MDAnalysis as mda
import numpy as np
import configparser as cp

import raijin_help

# Parse user input and options
ap = argparse.ArgumentParser(description=raijin_help.header,
                             formatter_class=argparse.RawDescriptionHelpFormatter,
                             usage="python Raijin.py [options]",
							 add_help=False,
							 exit_on_error=False)
ap._optionals.title = 'Options'
ap.add_argument('-h', '--help', action="store_true", help='Show this message and exit')
ap.add_argument('-top', type=str, default=None, required=False, metavar='',
                help='Topology file (.psf)')
ap.add_argument('-traj', type=str, default=None, required=False, metavar='',
                help='Trajectory file (.dcd)')
ap.add_argument('-outdir', type=str, default="Raijin_results", required=False, metavar='',
                help='Output folder in which results will be written (default: Raijin_results/)')
ap.add_argument('-config', type=str, default=None, required=False, metavar='',
                help='Input configuration file (default: config.dat)')
ap.add_argument('-template', type=str, default=None, required=False, metavar='',
                help='Create a template for input configuration file (default: config_template.dat)')

cmd = ap.parse_args()
start = timeit.default_timer()

if cmd.help is True:
	ap.print_help()
	print(raijin_help.help)
	sys.exit()
else:
	pass

top_file      = cmd.top
traj_file     = cmd.traj
outdir        = cmd.outdir
configinput   = cmd.config
template      = cmd.template
###############################################################################
# Stablishing default config values
config = cp.ConfigParser(allow_no_value=True, inline_comment_prefixes="#")

if template != None:
	with open(template, 'w') as configtemplate:
		configtemplate.write(raijin_help.template_content)
	print("\n>>> Configuration template file (" + str(template) + ") written.\n>>> Exiting...\n")
	sys.exit()
###############################################################################
# Reading configuration file and assuming absences
config.read(configinput)

try:
	sele_elecfield = str(config['Elecfield Selection']["sele_elecfield"].strip('"'))
except Exception as e1:
	sys.exit("\n>>> ERROR: sele_elecfield must be defined and a valid MDanalysis selection!\n")

try:
	mode = str(config['Target Selection']["mode"].strip('"')).lower()
	if mode == "atom":
		try:
			selatom = str(config['Target Selection']["selatom"].strip('"'))
		except Exception as e2:
			sys.exit("""\n>>> ERROR: in "atom" mode, selatom must be defined!\n""")
	elif mode == "bond":
		try:
			selbond1 = str(config['Target Selection']["selbond1"].strip('"'))
			selbond2 = str(config['Target Selection']["selbond2"].strip('"'))
		except Exception as e2:
			sys.exit("""\n>>> ERROR: in "bond" mode, both selbond1 and selbond2 must be defined!\n""")
	elif mode == "coordinate":
		try:
			tmp_targetcoordinate = str(config['Target Selection']["targetcoordinate"]).strip('[]').split(",")
			targetcoordinate = [float(item) for item in tmp_targetcoordinate]
		except Exception as e2:
			sys.exit("""\n>>> ERROR: in "coordinate" mode, a list of coordinates [X,Y,Z] must be provided!\n""")
	else:
		sys.exit("""\n>>> ERROR: "mode" must be defined as "atom", "bond" or coordinate!\n""")

except Exception as e2:
	sys.exit("""\n>>> ERROR: "mode" must be defined as "atom", "bond" or coordinate!\n""")


if mode == "coordinate":
	try:
		remove_self = str(config['Target Selection']["remove_self"]).lower()
		if remove_self  == "true":
			remove_self = True
			try:
				remove_cutoff = float(config['Target Selection']["remove_cutoff"])
			except Exception as e3:
				sys.exit("""\n>>> ERROR: if "remove_self" is True, a removal cutoff (float) is expected!\n""")
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
	else:
		include_solvent   = False
except:
	include_solvent = False

if include_solvent  == True:
	try:
		solvent_selection = str(config['Solvent']["solvent_selection"].strip('"'))
	except Exception as e4:
		sys.exit("""\n>>> ERROR: solvent_selection must be a valid MDanalysis selection!\n""")
	try:
		solvent_cutoff    = float(config['Solvent']["solvent_cutoff"])
	except Exception as e5:
		sys.exit("""\n>>> ERROR: solvent_cutoff must be a number!\n""")
else:
	pass

try:
	begintime = str(config['Time']["begintime"]).lower()
	if  begintime == "none":
		begintime = None
	else:
		begintime = int(begintime)
except:
	begintime = None

try:
	endtime = str(config['Time']["endtime"]).lower()
	if  endtime == "none":
		endtime = None
	else:
		endtime = int(endtime)
except:
	endtime = None

try:
	dt = int(config['Time']["dt"])
except:
	dt = 1


###############################################################################
# Being verbose about parameters chosen
print(raijin_help.header)
print("########################################################")
print(">>> Parameters used to run Raijin:")

print('[Elecfield Selection]')
print('sele_elecfield     = {}'.format(sele_elecfield))
print()
print("[Target Selection]")
print('mode               = {}'.format(mode))
if mode == "atom":
	print('selatom            = {}'.format(selatom))
elif mode == "bond":
	print('selbond1           = {}'.format(selbond1))
	print('selbond2           = {}'.format(selbond2))
elif mode == "coordinate":
	print('targetcoordinate   = {}'.format(targetcoordinate))
print('remove_self        = {}'.format(remove_self))
print()
print("[Solvent]")
print('include_solvent    = {}'.format(include_solvent))
if include_solvent == True:
	print('solvent_selection  = {}'.format(solvent_selection))
	print('solvent_cutoff     = {}'.format(solvent_cutoff))
print()
print("[Time]")
print('begintime          = {}'.format(begintime))
print('endtime            = {}'.format(endtime))
print('dt                 = {}'.format(dt))

###############################################################################
def mag(vector):
	""" Returns the magnitude of the vector.  """
	mag = np.linalg.norm(vector)
	return mag

def unit_vector(vector):
	""" Returns the unit vector of the vector.  """
	return vector / np.linalg.norm(vector)

def angle_between(v1, v2):

	v1_u = unit_vector(v1)
	v2_u = unit_vector(v2)

	angle_rad = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
	return angle_rad

def calc_EletricProperties(atom,refposition):

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

	if mode == "bond":
		#projection_u_on_v = (np.dot(u, v)/np.dot(v, v))*v
		proj_Ef_on_rbond = rbond_vec*(np.dot(Ef, rbond_vec)/np.dot(rbond_vec,rbond_vec))

		return Ef, proj_Ef_on_rbond
	else:
		pass

	return Ef

###############################################################################
# Opening output files
outdir = outdir + "/"
if os.path.exists(outdir):
	pass
else:
	os.system("mkdir " + outdir)

out = open(outdir + "ElecField.dat",'w')
out.write("""@    title "Electric Field"\n@    xaxis  label "Time (ps)"\n@    yaxis  label "MV/cm"\n""")
out.write("#time  Magnitude   Efield_X   Efield_Y   Efield_Z\n")
out.write("@type xy\n")

outres = open(outdir + "ElecField_per_residue.dat",'w')
outres.write("""@    title "Magnitude of Electric Field"\n@    xaxis  label "Residues"\n@    yaxis  label "MV/cm"\n""")
outres.write("#time  EField_Magnitude\n")
outres.write("@type xydy\n")

if mode == "bond":
	outproj = open(outdir + "ElecField_proj_onto_bond.dat",'w')
	outproj.write("""@    title "Electric Field Projection"\n@    xaxis  label "Time (ps)"\n@    yaxis  label "MV/cm"\n""")
	outproj.write("#time  Xcomp   Ycomp   Zcomp   Magnitude   Direction\n")
	outproj.write("@type xy\n")

	outalig = open(outdir + "ElecField_alignement.dat",'w')
	outalig.write("""@    title "Electric Field Projection"\n@    xaxis  label "Time (ps)"\n@    yaxis  label "Percentage"\n""")
	outalig.write("#time  Alignment_percentage\n")
	outalig.write("@type xy\n")

	outrbond = open(outdir + "bond_axis_info.dat", "w")
	outrbond.write("#time    rbondvec    rbondmag    rbondhat\n")

###############################################################################
# Creating Universe and making selections
if top_file != None and traj_file != None:
	u = mda.Universe(top_file, traj_file)
else:
	sys.exit("""\n>>> ERROR: Topology or Trajectory files missing. Are you not forgetting something?!\n""")

print("\n########################################################")
if mode == "atom":
	target_selection = u.select_atoms(selatom)
	if len(target_selection) > 1:
		print("\n>>> Running in ATOM mode!")
		print(">>> Target atoms (COG) = "+ str(target_selection.atoms) + "\n")
	else:
		print("\n>>> Running in ATOM mode!")
		print(">>> Target atom = "+ str(target_selection.atoms[0].name) +"-"+ str(target_selection.resnames[0]) + str(target_selection.resnums[0]) + "\n")
elif mode == "bond":
	selbond1 = u.select_atoms(selbond1)
	selbond2 = u.select_atoms(selbond2)
	target_selection = selbond1 + selbond2
	print("\n>>> Running in BOND mode!")
	print(">>> Target atoms = " + str(target_selection.atoms))
	print(">>> Bond axis direction = " + str(selbond1.atoms[0].name) +"-"+ str(selbond1.resnames[0]) + str(selbond1.resnums[0]) + " --> " + str(selbond2.atoms[0].name) +"-"+ str(selbond2.resnames[0]) + str(selbond2.resnums[0]) + "\n")
elif mode == "coordinate":
	target_selection = targetcoordinate
	print("\n>>> Running in COORDINATE mode!")
	print(">>> Target XYZ position = " +  str(target_selection))
	if remove_self == True:
		print(">>> Removing self contribution within a radial cutoff of " +  str(remove_cutoff) + " Angstroms\n")

else:
	sys.exit("\n>>> MODE should be 'bond' or 'atom' or 'coordinate'!")

# Making selection of atoms that exert the electric field
elecfield_selection = u.select_atoms(sele_elecfield)


###############################################################################
# Sanity check: atoms in target_selection should NOT be in elecfield_selection:
if mode == "atom" or mode == "bond":
	for atom in target_selection.atoms:
		if atom in elecfield_selection.atoms:
			sys.exit(">>> ERROR: Target atom(s) within ElecField selection! Consider using 'remove_self = True'.\n>>> Exiting...\n")

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
		dict_res_total[residue] = [[],[],[],[]]

###############################################################################
print("\n########################################################")
print("\n>>> Calculating Electric Field at time:")
if begintime == None:
	begintime = 0
if endtime == None:
	endtime = len(u.trajectory)
for ts in u.trajectory[begintime: endtime + dt:]:

	########################################################
	# Update environment and target selections
	if include_solvent == None:
		enviroment_selection = elecfield_selection + u.select_atoms("(around " + str(solvent_cutoff) + " protein) and " + solvent_selection)
	else:
		enviroment_selection = elecfield_selection

	if mode == "atom":
		if len(target_selection) > 1:
			refposition = target_selection.center_of_geometry()
		else:
			refposition = target_selection.atoms[0].position

	elif mode == "bond":
		position1 = selbond1.atoms[0].position
		position2 = selbond2.atoms[0].position

		rbond_vec = (position2 - position1) # axis from ref1 to ref2
		rbond_vec = rbond_vec*(10**(-10)) #convert from Angstrom to meter
		rbond_hat = unit_vector(rbond_vec)

		refposition = (position1 + position2)/2 #midway for both atoms

	elif mode == "coordinate":
		refposition = target_selection
	########################################################

	# Converting MDanalysis frames using defined dt
	frame = int(ts.frame) + 1
	time = "{:.0f}".format(frame*dt)
	print("Time = " + str(time) + " ps... (Target position: " + str(refposition) + ")")


	# Evaluate self_contribution removal
	if mode == "coordinate":
		# Take absolute coordinates from targetcoordinate
		coordX, coordY, coordZ = refposition
		# selects all atoms within a cutoff of a point in space (point X Y Z cutoff)
		self_contribution = u.select_atoms("point " + str(refposition[0]) + " " + str(refposition[1]) + " " + str(refposition[2]) + " " + str(remove_cutoff))

		if len(self_contribution) > 0:
			if remove_self == True:
				print(""">>> Warning: removing self contribution of: """ + str(self_contribution.atoms))
				# Remove self_contribution for this frame
				enviroment_selection = enviroment_selection - self_contribution
			else:
				print(""">>> Warning: some atoms are closed than """ + str(remove_cutoff) + """ A : """ + str(self_contribution.atoms))


	########################################################
	# opening a temporary dictionary to hold the contribution of each residue
	# for each frame
	dict_res_tmp = {}

	xfield_list = []
	yfield_list = []
	zfield_list = []
	xproj_list = []
	yproj_list = []
	zproj_list = []

	# Iterating over all atoms in selection to get their contributions
	for atom in enviroment_selection.atoms:
		residue = str(atom.resname) + "_" + str(atom.resid)
		if residue not in dict_res_tmp.keys() and residue in dict_res_total.keys():
			dict_res_tmp[residue] = [[],[],[]]
		else:
			pass

		if mode == "bond":
			Ef_xyz, proj_Ef_on_rbond = calc_EletricProperties(atom, refposition)
			xproj, yproj, zproj = proj_Ef_on_rbond
			xproj_list.append(xproj)
			yproj_list.append(yproj)
			zproj_list.append(zproj)
		else:
			Ef_xyz = calc_EletricProperties(atom, refposition)
		Efx, Efy, Efz = Ef_xyz

		xfield_list.append(Efx)
		yfield_list.append(Efy)
		zfield_list.append(Efz)

		# upload temp dict keys (residues) to include the contribution of each
		# atom in the residue. We will update the dictionary for each frame
		# later (line 470)
		if residue in dict_res_total.keys():
			res_tmp_x, res_tmp_y, res_tmp_z = dict_res_tmp[residue]
			res_tmp_x.append(Efx)
			res_tmp_y.append(Efy)
			res_tmp_z.append(Efz)
			dict_res_tmp[residue] = [res_tmp_x, res_tmp_y, res_tmp_z]
		else:
			pass

	# Sum the contributions of all atoms to get the resultant Efield
	totalEfx = np.sum(xfield_list)
	totalEfy = np.sum(yfield_list)
	totalEfz = np.sum(zfield_list)
	totalEf  = np.array([totalEfx, totalEfy, totalEfz])
	totalEfmag = mag(totalEf)

	# write the Efield output
	lineEfield  = str(time).ljust(10,' ') + str("{:.12e}".format(totalEfmag)).ljust(30,' ') + str("{:.12e}".format(totalEfx)).ljust(30,' ') + str("{:.12e}".format(totalEfy)).ljust(30,' ') + str("{:.12e}".format(totalEfz)).ljust(30,' ') + "\n"
	out.write(lineEfield)

	########################################################
	# Write information exclusive to BOND mode
	if mode == "bond":
		# Calculate projection
		totalprojx = np.sum(xproj_list)
		totalprojy = np.sum(yproj_list)
		totalprojz = np.sum(zproj_list)
		totalproj = np.array([totalprojx, totalprojy, totalprojz])
		# Calculate projection direction
		proj_direction = np.cos(angle_between(totalEf,rbond_hat))/abs(np.cos(angle_between(totalEf,rbond_hat)))
		totalproj = totalproj*proj_direction
		totalprojmag = mag(totalproj)
		# Write projection
		lineproj = str(time).ljust(10,' ') + str("{:.12e}".format(totalprojx)).ljust(30,' ') + str("{:.12e}".format(totalprojy)).ljust(30,' ') + str("{:.12e}".format(totalprojz)).ljust(30,' ') + str("{:.12e}".format(totalprojmag)).ljust(30,' ') + str(proj_direction) + "\n"
		outproj.write(lineproj)
		# Calculate and write percentage of projection alignment
		aligned = totalprojmag/totalEfmag
		linealig = str(time).ljust(10,' ') + str("{:.12e}".format(aligned)).ljust(30,' ') + "\n"
		outalig.write(linealig)

		# Write rbond_info data
		rvx, rvy, rvz = rbond_vec
		rhx, rhy, rhz = rbond_hat
		list1 = "[" + str(rvx) + "," + str(rvy) + "," + str(rvz) + "]"
		list2 = "[" + str(rhx) + "," + str(rhy) + "," + str(rhz) + "]"
		outrbond.write(str(time).ljust(10,' ') +  list1  + "        " + str(mag(rbond_vec)) + "        " + list2 + "\n")

	########################################################
	# Update dictionary with the contribution of each residue with the
	# values of Efield of each residue in this frame
	for r, comp in dict_res_tmp.items():
		res_efield_x, res_efield_y, res_efield_z = comp

		resEfx = np.sum(res_efield_x)
		resEfy = np.sum(res_efield_y)
		resEfz = np.sum(res_efield_z)
		resEfmag = mag(np.array([resEfx, resEfy, resEfz]))

		resEfx_tmp, resEfy_tmp, resEfz_tmp, resEfmag_tmp = dict_res_total[r]
		resEfx_tmp.append(resEfx)
		resEfy_tmp.append(resEfy)
		resEfz_tmp.append(resEfz)
		resEfmag_tmp.append(resEfmag)
		dict_res_total[r] = [resEfx_tmp, resEfy_tmp, resEfz_tmp, resEfmag_tmp]

###############################################################################
# Calculate the average contribution of each residue throughout trajectory
for r, comp in dict_res_total.items():
	r = r.partition('_')[2] # ignore resname and keep resid
	res_efield_x, res_efield_y, res_efield_z, resEfmag = comp
	resEfmag_avg = np.average(resEfmag)
	resEfmag_std = np.std(resEfmag)

	# write contribution of each residue
	lineres = str(r).ljust(10,' ') + str("{:.12e}".format(resEfmag_avg)).ljust(30,' ') + str("{:.12e}".format(resEfmag_std)).ljust(30,' ') + "\n"
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
