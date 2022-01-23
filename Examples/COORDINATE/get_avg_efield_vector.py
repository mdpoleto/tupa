import os
import numpy as np


####################################
def make_avg(file):

	f = open(file, 'r')

	Xlist = []
	Ylist = []
	Zlist = []

	while True:
		line = f.readline()
		if not line: break

		if line[0] != "#" and line[0] != "@":
			Xvalue = float(line.split()[2])
			Yvalue = float(line.split()[3])
			Zvalue = float(line.split()[4])
			Xlist.append(Xvalue)
			Ylist.append(Yvalue)
			Zlist.append(Zvalue)

	avgX = np.average(Xlist)
	avgY = np.average(Ylist)
	avgZ = np.average(Zlist)

	return avgX, avgY, avgZ


######
for i in range(26,73):

	i = str(i)
	file = "Z_" + i + "/ElecField.dat"

	avgX, avgY, avgZ = make_avg(file)

	line = "efield_point [0,0," + i +"], efield=[" + str(avgX) + ", " + str(avgY) + ", " + str(avgZ) + "], scale=0.01, radius=0.3, color=teal"
	print(line)

######
