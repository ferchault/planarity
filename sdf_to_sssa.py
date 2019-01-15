#!/usr/bin/env python
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import scipy.spatial

def get_sssa_components(coordinates):
    hull = scipy.spatial.ConvexHull(coordinates, qhull_options='QJ')   
    #diameter = max(scipy.spatial.distance.pdist(coordinates[hull.vertices], 'euclidean'))

    return hull.volume, hull.area

fn = sys.argv[1]
for mol in Chem.SDMolSupplier(fn):
	n = mol.GetNumAtoms()
	if AllChem.EmbedMolecule(mol) == -1:
		continue
	V, A = get_sssa_components(mol.GetConformer().GetPositions())
	print (V, A, n)
