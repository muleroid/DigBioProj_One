# Rewrite Shen's parse_pdb.py using Prody

import prody as pr

#______________________________________________________________
# run through a single pdb file
#def run(pdbFile):
def getHforAtom(aParsedProPDB, anAtom):
	for index in range(anAtom.getIndex()-5, anAtom.getIndex()+5):
		otherAtom = aParsedPDBfile[index]
		# 1-BE HYDROGEN, 2-BE SAME RESNUM, 3-
		if(otherAtom.getElement() == 'H'):
			if(otherAtom.getResnum() == anAtom.getResnum()):
				#THIS GIVES MANY
				a=1 
	#return pr.findNeighbors(anAtom,1)
	#hydrogens = aParsedPDBfile.select('element H')
	#hydrogens.select('')

pdbFile = '/Users/fsimon/Google Drive/School/Winter_13/Digital Biology/Project/DigBioProj_One/pdbs/1A23_A.pdb'

aParsedPDBfile = pr.parsePDB(pdbFile, model=1, secondary=True, chain='A', altLoc=False)

# get secondary structure by aParsedPDBfile[i].getSecstr()

nitrogens = aParsedPDBfile.select('element N')
oxygens = aParsedPDBfile.select('element O')
hydrogens = aParsedPDBfile.select('element H')

# antecedents = 
donors = nitrogens
acceptors = oxygens
'''for donor in donors:
	for acceptor in acceptors:
		dist = pr.calcDistance(donor, acceptor)
		if(dist >= 3.5):
			continue
		hydrogen = getHforAtom(aParsedPDBfile, acceptor)
		hdistance = pr.calcDistance(acceptor, hydrogen)
		if(hdistance >= 2.5):
			continue
		dha = pr.calcAngle(donor, hydrogen, acceptor)
		if(dha <= 90):
			conitnue
		#We have a valid h-bond
		print donor, acceptor'''
