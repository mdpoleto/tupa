import os, sys
import numpy as np

#file = sys.argv[1]

####################################
def make_avg(file):

        f = open(file, 'r')

        maglist = []
        Xlist = []
        Ylist = []
        Zlist = []

        while True:
                line = f.readline()
                if not line: break

                if line[0] != "#" and line[0] != "@":
                        mag    = float(line.split()[1])
                        Xvalue = float(line.split()[2])
                        Yvalue = float(line.split()[3])
                        Zvalue = float(line.split()[4])
                        maglist.append(mag)
                        Xlist.append(Xvalue)
                        Ylist.append(Yvalue)
                        Zlist.append(Zvalue)

        mag  = np.average(maglist)
        std  = np.std(maglist)
        avgX = np.average(Xlist)
        avgY = np.average(Ylist)
        avgZ = np.average(Zlist)

        return mag, std, avgX, avgY, avgZ


######

for i in range(2,32,2):

    mag, std, avgX, avgY, avgZ = make_avg("bla_" + str(i) + "/ElecField.dat")

    line = "efield_point ATOM, efield=[" + str(avgX) + ", " + str(avgY) + ", " + str(avgZ) + "], scale=0.01, radius=0.3, color=teal"
    #print(line)
    print(str(i) + "    " + str(mag) + "    " + str(std))

######
