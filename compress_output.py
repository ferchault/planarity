#!/usr/bin/env python
""" Reduces a file with detailed output (ID, volume, area, distance, numatoms) into a binary numpy file with S and area/atom as only columns."""
import sys
import numpy as np

infile, outfile = sys.argv[1:]
fh = open(infile)
collected = []
for line in fh:
	parts = line.strip().split()
	if len(parts) != 5:
		continue
	V, A, d, N = map(float, parts[1:])
	collected.append((A*d/V, A/N))

np.save(outfile, np.array(collected))
