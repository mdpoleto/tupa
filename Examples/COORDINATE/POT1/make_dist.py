import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from scipy import stats
import sys, os


file = sys.argv[1]
outname = sys.argv[2]

##################################################

def make_distribution(given_list):
	dist = stats.gaussian_kde(given_list)

	minvalue = min(given_list)*1.5# - min(given_list)*0.2
	maxvalue = max(given_list)*1.5# + max(given_list)*0.2


	ord_values = np.linspace(-200,200,1000)
	dist = list(dist(ord_values))

	raw_list = []

	for i in range(len(ord_values)):
		dip = ord_values[i]
		pdf = dist[i]
		raw_list.append([dip,pdf])

	return raw_list

##################################################
f = open(file, 'r')

value_list = []


print(">>> Parsing elements...")
line = f.readline() #skip the header
while True:
	line = f.readline()
	if not line: break

	if line[0] == "#" or line[0] == "!" or line[0] == "@":
		pass
	else:
		frame = line.split()[0]
		value = line.split()[4]
		value_list.append(float(value))



outdist = open(outname, "w")

print(">>> Calculating distribution...")
dist = make_distribution(value_list)

for d in dist:
	val = str(d[0])
	pdf = str(d[1])
	outdist.write(val.ljust(30, ' ') + pdf.ljust(10,' ') + "\n")


outdist.close()
print(">>> Finished!")
