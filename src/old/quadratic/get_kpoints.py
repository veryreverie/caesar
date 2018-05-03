#! /usr/bin/env python
import sys
fileName = sys.argv[1]
lines = open(fileName).read().splitlines()

# Get k-points
lookup = 'number of k points='
line_with_lookup = [l for l in lines if lookup in l][0]
index_line_with_lookup = lines.index(line_with_lookup)

number_of_k_points = int(line_with_lookup.split()[4])

#print "Reading file... "
newFile_kpoints = [];
for iline in range(index_line_with_lookup+2, index_line_with_lookup+number_of_k_points+2):
	thisLine = lines[iline].split()
	k1 = thisLine[-6]
	k2 = thisLine[-5]
	k3 = thisLine[-4][0:-2]
	wk = thisLine[-1]
	newLine = "{0} {1} {2} {3} ".format(k1, k2, k3, wk)
	newFile_kpoints.append(newLine)
#print "Done."


# Get bands
newFile_bands = []
first_line = [l for l in lines if " k =" in l][0]
idx_first_line = lines.index(first_line)+2
idx_lines_with_bands = range(0, number_of_k_points)
idx_lines_with_bands = [x*7 for x in idx_lines_with_bands]
idx_lines_with_bands = [x+idx_first_line for x in idx_lines_with_bands]
newFile_bands = [lines[i] for i in idx_lines_with_bands]


# Save file
if len(newFile_kpoints)!=len(newFile_bands):
	sys.exit("Oops! Something happened.")

newFile = open('kpt_band.dat','write')
for i in range(0, number_of_k_points):
	newLine = "{0} {1}\n".format(newFile_kpoints[i], newFile_bands[i])
	newFile.write(newLine)
newFile.close()
#print "Saved file kpt.dat"
