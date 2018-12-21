#!/usr/bin/env python
import scipy.spatial
import numpy as np
import sys

def get_sssa_components(coordinates):
    hull = scipy.spatial.ConvexHull(coordinates, qhull_options='QJ')   
    diameter = max(scipy.spatial.distance.pdist(coordinates[hull.vertices], 'euclidean'))

    return {'diameter': diameter, 'volume': hull.volume, 'area': hull.area}

fn, identifier = sys.argv[1:]
fh = open(fn)
molline = 0
coordinates = []
name = 'noname'
for line in fh:
	molline += 1
	if molline == 4:
		numatoms = int(line[:3].strip())
	if molline > 4 and numatoms > 0:
		coordinates.append(list(map(float, line.strip().split()[:3])))
		numatoms -= 1
	if molline > 4 and numatoms == 0:
		if name is None:
			name = line.strip()
		if line.strip() == '> <%s>' % identifier:
			name = None
		if line.strip() == '$$$$':
			res = get_sssa_components(np.array(coordinates))
			print (name, res['volume'], res['area'], res['diameter'])
			coordinates = []
			molline = 0
			name = 'noname'
