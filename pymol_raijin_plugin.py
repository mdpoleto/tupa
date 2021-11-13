from chempy import cpv
from pymol import cmd, cgo, CmdException
import numpy as np

'''
CREATION
	(c) Marcelo D. Poleto, Sep 2021. Virginia Tech
	This script was based on 'cgo_arrow.py' by Thomas Holder.
DESCRIPTION
	Allows the user to create arrows representing:
	1) bond dipoles between 2 selected atoms (atom1 -> atom2)
	2) Electric field vectors midway between 2 picked atoms ([atom1+atom2]/2)
	3) Electric field vector at a given atom or coordinate
ARGUMENTS
	bond_atom1 = string: single atom selection or list of 3 floats {default: pk1}
	bond_atom2 = string: single atom selection or list of 3 floats {default: pk2}
	efield = list of 3 floats containing the XYZ electric field components {default: [1.0. 1.0, 1.0]}
	radius = float: arrow radius {default: 0.1}
	scale = float: scale factor to change arrow size {default: 0.0}
	hlength = float: length of arrow head in percentage of efield magnitude {default: 30%}
	hradius = float: radius of arrow head in percentage of radius {default: 200%}
	color = string: one or two color names {default: blue red}
	name = string: name of CGO object
'''

def mag(vector):
	mag = np.linalg.norm(vector)
	return mag

def get_coord(v):
	if not isinstance(v, str):
		return v
	if v.startswith('['):
		return cmd.safe_list_eval(v)
	return cmd.get_atom_coords(v)

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

	print("\n##############################")
	print("Bond dipole unit vectors (r_hat)= ", normal)
	print("##############################\n")

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

def efield_bond(bond_atom1='pk1', bond_atom2='pk2', efield=[1.0, 1.0, 1.0], scale=1.0, radius=0.1, hlength=0.3, hradius=None, color='blue red', name=''):

	radius, scale, hlength = float(radius), float(scale), float(hlength)
	if not 0 <= hlength <= 1.0:
		print("\nERROR! hlength must be between 0 and 1!")
		return
	else:
		hlength = 1 - hlength

	if hradius == None:
		hradius = float(radius)*2 # arrow head 2 times thicker
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
	xyz = cpv.average(xyz1, xyz2)

	def get_vec(v):
		if v.startswith('['):
			v = v.strip("[]").split(",")
			vecx = float(v[0])
			vecy = float(v[1])
			vecz = float(v[2])
			vector = [vecx, vecy, vecz]
			return vector

	efield = get_vec(efield)
	efield = cpv.scale(efield,scale)
	efieldhat = cpv.normalize(efield) # create the unit vector of vector efield

	print("\n##############################")
	print("Electric Field (scaled)= ", efield)
	print()
	print("Electric Field unit vectors (Ef_hat)= ", efieldhat)
	print("##############################\n")

	# The arrow is a cylinder plus a cone. So we actually draw a cylinder and
	# discount from it the size of the arrow head (the cone).
	end = cpv.add(efield,xyz)
	xyz_cyl = cpv.add(cpv.scale(efield, hlength), xyz)

	obj = [cgo.CYLINDER] + xyz + xyz_cyl + [radius] + color1 + color2 + \
		[cgo.CONE] + xyz_cyl + end + [hradius, 0.0] + color2 + color2 + \
		[1.0, 0.0]

	if not name:
		name = cmd.get_unused_name('efield')

	cmd.load_cgo(obj, name)

def efield_point(point='pk1', efield=[1.0, 1.0, 1.0], scale=1.0, radius=0.1, hlength=0.3, hradius=None, color='blue red', name=''):

	radius, scale, hlength = float(radius), float(scale), float(hlength)
	if not 0 <= hlength <= 1.0:
		print("\nERROR! hlength must be between 0 and 1!")
		return
	else:
		hlength = 1 - hlength

	if hradius == None:
		hradius = float(radius)*2 # arrow head 2 times thicker
	else:
		hradius = float(hradius)

	try:
		color1, color2 = color.split()
	except:
		color1 = color2 = color
	color1 = list(cmd.get_color_tuple(color1))
	color2 = list(cmd.get_color_tuple(color2))

	xyz = get_coord(point)

	def get_vec(v):
		if v.startswith('['):
			v = v.strip("[]").split(",")
			vecx = float(v[0])
			vecy = float(v[1])
			vecz = float(v[2])
			vector = [vecx, vecy, vecz]
			return vector

	efield = get_vec(efield)
	efield = cpv.scale(efield,scale)
	efieldhat = cpv.normalize(efield) # create the unit vector of vector efield

	print("\n##############################")
	print("Electric Field (scaled)= ", efield)
	print()
	print("Electric Field unit vectors (Ef_hat)= ", efieldhat)
	print("##############################\n")

	# The arrow is a cylinder plus a cone. So we actually draw a cylinder and
	# discount from it the size of the arrow head (the cone).
	end = cpv.add(efield,xyz)
	xyz_cyl = cpv.add(cpv.scale(efield, hlength), xyz)

	obj = [cgo.CYLINDER] + xyz + xyz_cyl + [radius] + color1 + color2 + \
		[cgo.CONE] + xyz_cyl + end + [hradius, 0.0] + color2 + color2 + \
		[1.0, 0.0]

	if not name:
		name = cmd.get_unused_name('efield')

	cmd.load_cgo(obj, name)
###################################################

cmd.extend('draw_bond_axis', draw_bond_axis)
cmd.extend('efield_bond', efield_bond)
cmd.extend('efield_point', efield_point)
