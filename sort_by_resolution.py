#Sort pdbs into seperate resolution directories

import os, re, shutil, StringIO

directory = "/Users/fsimon/Desktop/nrpdb/unbucketed/"
dirbuck0  = "/Users/fsimon/Desktop/nrpdb/buck0/"
dirbuck1  = "/Users/fsimon/Desktop/nrpdb/buck1/"
dirbuck2  = "/Users/fsimon/Desktop/nrpdb/buck2/"
dirbuck3  = "/Users/fsimon/Desktop/nrpdb/buck3/"
dirbuck4  = "/Users/fsimon/Desktop/nrpdb/buck4/"
dirbuck5  = "/Users/fsimon/Desktop/nrpdb/buck5/"
dirbuck6  = "/Users/fsimon/Desktop/nrpdb/buck6/"
dirbuck7  = "/Users/fsimon/Desktop/nrpdb/buck7/"
dirbuck8  = "/Users/fsimon/Desktop/nrpdb/buck8/"
dirbuck_nores  = "/Users/fsimon/Desktop/nrpdb/buck_nores/"


allresolutions = []

int1 = 0
for item in os.listdir(directory):
	int1 = int1+1
	print int1
	resolution = -1
	if(item.startswith(".")):
		continue
	if(not item.endswith("pdb")):
		continue
	fil = open(directory+item)
	file_lines = fil.readlines()
	for line in file_lines:
		if ("2 RESOLUTION. NOT APPLICABLE." in line):
			#print resolution
			continue
		m = re.search("RESOLUTION. ([0-9]*.[0-9]+) ANGSTROMS.", line)
		if m:
			resolution = float(m.groups()[0])
			#print resolution
			continue
	allresolutions.append(float(resolution))
	#Now we have the resolution.
	if(resolution == -1):
		wrfi = open(dirbuck_nores+item, 'w')
		wrfi.writelines(file_lines)
		wrfi.close()
		fil.close()
	if((resolution > 0.5)) & ((resolution < 1.0)):
		wrfi = open(dirbuck0+item, 'w')
		wrfi.writelines(file_lines)
		wrfi.close()
		fil.close()
	if((resolution > 1)) & ((resolution < 1.25)):
		wrfi = open(dirbuck1+item, 'w')
		wrfi.writelines(file_lines)
		wrfi.close()
		fil.close()	
	if((resolution > 1.25) & (resolution < 1.5)):
		wrfi = open(dirbuck2+item, 'w')
		wrfi.writelines(file_lines)
		wrfi.close()
		fil.close()	
	if((resolution > 1.5) & (resolution < 1.75)):
		wrfi = open(dirbuck3+item, 'w')
		wrfi.writelines(file_lines)
		wrfi.close()
		fil.close()	
	if((resolution > 1.75) & (resolution < 2)):
		wrfi = open(dirbuck4+item, 'w')
		wrfi.writelines(file_lines)
		wrfi.close()
		fil.close()	
	if((resolution > 2) & (resolution < 2.25)):
		wrfi = open(dirbuck5+item, 'w')
		wrfi.writelines(file_lines)
		wrfi.close()
		fil.close()
	if((resolution > 2.25) & (resolution < 2.5)):
		wrfi = open(dirbuck6+item, 'w')
		wrfi.writelines(file_lines)
		wrfi.close()
		fil.close()	
	if((resolution > 2.5) & (resolution < 2.75)):
		wrfi = open(dirbuck7+item, 'w')
		wrfi.writelines(file_lines)
		wrfi.close()
		fil.close()	
	if((resolution > 2.75) & (resolution < 3)):
		wrfi = open(dirbuck8+item, 'w')
		wrfi.writelines(file_lines)
		wrfi.close()
		fil.close()	
	if(resolution > 3):
		print "GREATER THAN 3"