import os, subprocess

#Run MolProbity on Examples
path = '/Users/fsimon/Google Drive/School/Winter_13/Digital Biology/Project/DigBioProj_One/pdbs/'
pdbs = os.listdir(path)
for pdb in pdbs:
	#Add Hydrogens with REDUCE
	input = ["/Users/fsimon/Downloads/reduce", "-build", path+pdb, '>', path+pdb.replace('.pdb','')+'H.pdb']
	print input
	a = subprocess.call(input)
	print a