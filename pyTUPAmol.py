#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#Marcelo D. Poleto
#FEB 2022

from chempy import cpv
from pymol import cmd, cgo, CmdException
import numpy as np
import math

'''
CREATION
	(c) Marcelo D. Poleto, Feb 2022. Virginia Tech
	This script was based on 'cgo_arrow.py' by Thomas Holder.
DESCRIPTION
	Allows the user to create arrows representing:
	1) bond dipoles between 2 selected atoms (atom1 -> atom2)
	2) Electric field vectors midway between 2 picked atoms ([atom1+atom2]/2)
	3) Electric field vector at a given atom or coordinate
ARGUMENTS
	bond_atom1 = string: single atom selection or list of 3 floats {default: pk1}
	bond_atom2 = string: single atom selection or list of 3 floats {default: pk2}
	point      = string: single atom selection or list of 3 floats {default: pk1}

	efield      = list of 3 floats containing the XYZ electric field components {default: [1.0. 1.0, 1.0]}
	radius      = float: arrow radius {default: 0.1}
	scale       = float: scale factor to change arrow size {default: 0.0}
	hlength     = float: length of arrow head in percentage of efield magnitude {default: 30%}
	hradius     = float: radius of arrow head in percentage of radius {default: 2*radius}
	color       = string: one or two color names {default: blue red}
	stdev       = angle to define the spatial standard deviation of the efield
	efield_name = string: name of CGO object for the efield vector
	stdev_name  = string: name of CGO object for the 3D standard deviation
'''

def mag(vector):
	mag = np.linalg.norm(vector)
	return mag

def get_vec(array):
	if array.startswith('['):
		v = array.strip("[]").split(",")
		vecx = float(v[0])
		vecy = float(v[1])
		vecz = float(v[2])
		vector = [vecx, vecy, vecz]
		return vector
	else:
		print("Warning! Could not make a vector from the efield array!")
		return False

def get_coord(v):
	# An ugly workaround to issue an error in case neither cases are contemplated
	position = False

	try:
		# if it is a coordinate array
		if v.startswith('['):
			position = cmd.safe_list_eval(v)
			return position
		else:
			pass
	except:
		pass

	try:
		# if it is a pk1 or a syntax selection
		nAtoms = cmd.count_atoms(v)
		if nAtoms==1:
			# atom coordinates
			position = cmd.get_atom_coords(v)
			return position
		elif nAtoms > 1:
			# more than one atom --> use selection center of geometry
			centroid = cpv.get_null()
			model = cmd.get_model(v)
			for a in model.atom:
				centroid = cpv.add(centroid, a.coord)

			position = cpv.scale(centroid, 1. / nAtoms)
			return position
		else:
			pass
	except:
		pass

	return position


def draw_bond_axis(atom1='pk1', atom2='pk2', radius=0.1, gap=0.5, hlength=0.4, hradius=0.2, color='gray60', name=''):
	radius, gap = float(radius), float(gap)
	hlength, hradius = float(hlength), float(hradius)

	try:
		color1, color2 = color.split()
	except:
		color1 = color2 = color
	color1 = list(cmd.get_color_tuple(color1))
	color2 = list(cmd.get_color_tuple(color2))

	xyz1 = get_coord(atom1)
	xyz2 = get_coord(atom2)
	normal = cpv.normalize(cpv.sub(xyz1, xyz2))

	print("\n###############################")
	print("###### Running pyTUPÃmol ######")
	print("Bond axis unit vectors (r_hat)= ", normal)
	print("###############################\n")

	if hlength < 0 or hradius < 0:
		hlength=0.2
		hradius=0.2

	if gap:
		diff = cpv.scale(normal, gap)
		xyz1 = cpv.sub(xyz1, diff)
		xyz2 = cpv.add(xyz2, diff)

	xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

	obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
		[cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
		[1.0, 0.0]

	if not name:
		name = cmd.get_unused_name('bond_dipole')

	cmd.load_cgo(obj, name)
	return

def efield_bond(bond_atom1='pk1', bond_atom2='pk2', efield=None, scale=1.0, radius=0.1, hlength=0.3, hradius=None, color='blue red', stdev=0.0, efield_name='', stdev_name=''):

	if efield == None:
		efield = list([1.0, 1.0, 1.0])
	else:
		pass

	radius, scale, hlength, stdev = float(radius), float(scale), float(hlength), float(stdev)
	if not 0 <= hlength <= 1.0:
		print("\nERROR! hlength must be between 0 and 1!")
		return
	else:
		hlength = 1 - hlength

	if hradius == None:
		hradius = float(radius)*2 # arrow head 2 times more width than the cylinder
	else:
		hradius = float(hradius)

	try:
		color1, color2 = color.split()
	except:
		color1 = color2 = color
	color1 = list(cmd.get_color_tuple(color1))
	color2 = list(cmd.get_color_tuple(color2))

	xyz1 = get_coord(bond_atom1) # get the coordinates from which the arrow will be drawn
	xyz2 = get_coord(bond_atom2) # get the coordinates which the arrow will be drawn to

	if xyz1 == False or xyz2 == False:
		print("\nERROR! Could not find coordinates for the provided selection!")
		return
	else:
		xyz = cpv.average(xyz1, xyz2)

	efield = get_vec(efield)
	efield = cpv.scale(efield,scale)
	efieldhat = cpv.normalize(efield) # create the unit vector of vector efield
	efieldmag = mag(efield)

	print("\n##############################")
	print("##### Running pyTUPÃmol ######")
	print("Probe position = ", xyz)
	print("Electric Field (scaled) = ", efield)
	print("Electric Field magnitude = ", efieldmag)
	print("Electric Field unit vectors (Ef_hat) = ", efieldhat)
	print("##############################\n")

	# The arrow is a cylinder plus a cone. So we actually draw a cylinder and
	# discount from it the size of the arrow head (the cone).
	end = cpv.add(efield,xyz)
	xyz_cyl = cpv.add(cpv.scale(efield, hlength), xyz)

	obj = [cgo.CYLINDER] + xyz + xyz_cyl + [radius] + color1 + color2 + \
		[cgo.CONE] + xyz_cyl + end + [hradius, 0.0] + color2 + color2 + \
		[1.0, 0.0]

	if not efield_name:
		efield_name = cmd.get_unused_name('efield')
	cmd.load_cgo(obj, efield_name)


	if stdev > 0.0:
		stdev_rad = math.radians(stdev)
		std_radius = mag(cpv.scale(xyz_cyl,math.sin(stdev_rad)))
		std_height = cpv.scale(xyz_cyl,math.cos(stdev_rad))

		obj2 = [cgo.ALPHA, 0.25] + [cgo.CONE] + xyz + std_height + [radius, std_radius] + color1 + color2 + [1.0, 1.0]

		if not stdev_name:
			stdev_name = cmd.get_unused_name('spatialdev')

		cmd.load_cgo(obj2, stdev_name)

def efield_point(point='pk1', efield=None, scale=1.0, radius=0.1, hlength=0.3, hradius=None, color='blue red', stdev=0.0, efield_name='', stdev_name=''):

	if efield == None:
		efield = list([1.0, 1.0, 1.0])
	else:
		pass

	radius, scale, hlength, stdev = float(radius), float(scale), float(hlength), float(stdev)
	if not 0 <= hlength <= 1.0:
		print("\nERROR! hlength must be between 0 and 1!")
		return
	else:
		hlength = 1 - hlength

	if hradius == None:
		hradius = float(radius)*2 # arrow head 2 times more width than the cylinder
	else:
		hradius = float(hradius)

	try:
		color1, color2 = color.split()
	except:
		color1 = color2 = color
	color1 = list(cmd.get_color_tuple(color1))
	color2 = list(cmd.get_color_tuple(color2))

	xyz = get_coord(point)
	if xyz == False:
		print("\nERROR! Could not find coordinates for the provided selection!")
		return

	efield = get_vec(efield)
	efield = cpv.scale(efield,scale)
	efieldhat = cpv.normalize(efield) # create the unit vector of vector efield
	efieldmag = mag(efield)

	print("\n##############################")
	print("##### Running pyTUPÃmol ######")
	print("Probe position = ", xyz)
	print("Electric Field (scaled) = ", efield)
	print("Electric Field magnitude = ", efieldmag)
	print("Electric Field unit vectors (Ef_hat) = ", efieldhat)
	print("##############################\n")

	# The arrow is a cylinder plus a cone. So we actually draw a cylinder and
	# discount from it the size of the arrow head (the cone).
	end = cpv.add(efield,xyz)
	xyz_cyl = cpv.add(cpv.scale(efield, hlength), xyz)

	obj = [cgo.CYLINDER] + xyz + xyz_cyl + [radius] + color1 + color2 + \
		[cgo.CONE] + xyz_cyl + end + [hradius, 0.0] + color2 + color2 + \
		[1.0, 0.0]

	if not efield_name:
		efield_name = cmd.get_unused_name('efield')
	cmd.load_cgo(obj, efield_name)


	if stdev > 0.0:
		stdev_rad = math.radians(stdev)
		std_radius = mag(cpv.scale(xyz_cyl,math.sin(stdev_rad)))
		std_height = cpv.scale(xyz_cyl,math.cos(stdev_rad))

		obj2 = [cgo.ALPHA, 0.25] + [cgo.CONE] + xyz + std_height + [radius, std_radius] + color1 + color2 + [1.0, 1.0]

		if not stdev_name:
			stdev_name = cmd.get_unused_name('spatialdev')

		cmd.load_cgo(obj2, stdev_name)
###################################################

cmd.extend('draw_bond_axis', draw_bond_axis)
cmd.extend('efield_bond', efield_bond)
cmd.extend('efield_point', efield_point)
