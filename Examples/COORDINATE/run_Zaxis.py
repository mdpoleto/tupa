import os
import numpy as np


####################################
def make_avg(file):

	f = open(file, 'r')

	values_list = []

	while True:
		line = f.readline()
		if not line: break

		if line[0] != "#" and line[0] != "@":
			Zvalue = float(line.split()[4])
			values_list.append(Zvalue)

	avg = np.average(values_list)
	std = np.std(values_list)

	return avg,std

def make_config_file(filename, iteration):

	o = open(filename, 'w')

	content = """
[Elecfield Selection]
# The atoms from which we calculate the electric field
sele_elecfield      = segid PRO* or segid MEMB

[Target Selection]
# Provide the target selection for the MODE of you choice
# e.g. if bond is used, then modebond1 and modebond2 must be defined.
mode                = COORDINATE    # ATOM or BOND or COORDINATE
targetcoordinate    = [0,0,{iter}]
remove_self         = False   # For COORDINATE mode only, wether remove the
                              # contribution of self within a cutoff of the coordinate

[Solvent]
include_solvent     = False    # or False
solvent_cutoff      = 10
solvent_selection   = segid TIP3

[Time]
dt                  = 10      # Frequency of frames written in your trajectory (in picosecond)
""".format(iter = str(iteration))

	o.write(content)
	o.close()
####################################


######
o = open("efield_Zaxis.dat", 'w')

for i in range(26,73):

	i = str(i)
	make_config_file("config.conf", i)

	if not os.path.exists("Z_" + i):
		os.system("TUPA.py -top 3ouf_c36.equil.psf -traj 3ouf_c36.traj.dcd -conf config.conf -outdir Z_" + i)
######

######
for i in range(26,73):

	i = str(i)
	file = "Z_" + i + "/ElecField.dat"

	avg, std = make_avg(file)

	line = str(i).ljust(10, ' ') + str(avg).ljust(25, ' ') + str(std).ljust(10,' ') + "\n"

	o.write(line)

o.close()
######
