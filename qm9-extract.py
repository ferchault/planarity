#!/usr/bin/env python
import scipy.spatial
import numpy as np
import sys
import glob

def get_sssa_components(coordinates):
    hull = scipy.spatial.ConvexHull(coordinates, qhull_options='QJ')
    return hull.volume, hull.area

def coordinate_array(fn):
	lines = open(fn).readlines()
	numatoms = int(lines[0])
	coords = []
	for atom in range(numatoms):
		coords.append(list(map(float, lines[2+atom].strip().split()[1:4])))
	return np.array(coords), lines[1].strip()

fn = '/mnt/c/Users/guido/data/qm9/coord/46/3e/dsgdb9nsd_086034.xyz'
def do_fn(fn):
	c, label = coordinate_array(fn)
	print (label, *get_sssa_components(c))

for xyzfile in glob.glob('/mnt/c/Users/guido/data/qm9/coord/*/*/*.xyz'):
	try:
		do_fn(xyzfile)
	except:
		continue
