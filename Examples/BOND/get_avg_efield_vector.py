import os, sys
import numpy as np

file = sys.argv[1]

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

avgX, avgY, avgZ = make_avg(file)

line = "efield_bond segid HETA and name C1, segid HETA and name O1, efield=[" + str(avgX) + ", " + str(avgY) + ", " + str(avgZ) + "], scale=0.01, radius=0.3, color=cyan"
print(line)

######
