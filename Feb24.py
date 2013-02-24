import prody as pr
import numpy as np

DONOR_ACCEPTOR_MAXDISTANCE = 3.5
HYDROGEN_ACCEPTOR_MAXDISTANCE = 2.5
DHA_ANGLE = 90
DAB_ANGLE = 90
HAB_ANGLE = 90

#COLUMNS FOR DATA
#STRUCTURE INDICES ARE: 0=ALPHA, 1=310ALPHA, 2=PARALLEL, 3=ANTIPARALLEL
#COLUMN[STRUCTURE][1-N]
COLUMN_D_ON  = [[],[],[]] # Convert to [[],[],[],[]] when we separate parallel and antiparallel
COLUMN_D_OH  = [[],[],[]] # Convert to [[],[],[],[]] when we separate parallel and antiparallel
COLUMN_A_NHO = [[],[],[]] # Convert to [[],[],[],[]] when we separate parallel and antiparallel
COLUMN_A_HOC = [[],[],[]] # Convert to [[],[],[],[]] when we separate parallel and antiparallel
COLUMN_BETA  = [[],[],[]] # Convert to [[],[],[],[]] when we separate parallel and antiparallel
COLUMN_GAMMA = [[],[],[]] # Convert to [[],[],[],[]] when we separate parallel and antiparallel


def getHforAtom(aParsedProPDB, anAtom):
	resnumneighbors = aParsedProPDB.select('resnum '+str(anAtom.getResnum()))
	for at in resnumneighbors:
		if(at.getElement() == 'H'):
			#FILTER OTHER HYDROGENS
			return at

def getAntecedent (apdb, anAtom):
	aminoGroup = apdb.select('resnum ' + str(anAtom.getResnum()))
	for at in aminoGroup:
		if(at.getName() == 'C'):
			return at

def getSSIndex(acc_ante):
	#RETURN 0=ALPHA, 1=310ALPHA, 2=PARALLEL, 3=ANTIPARALLEL
	#DSSP (Pauli) structure convention
	ante_str = acc_ante.getSecstr()
	if( ante_str == 'H' ):
		#alpha
		return 0
	if( ante_str == 'G' ):
		#310
		return 1
	if( ante_str == 'E' ):
		#BOTH PARALLEL AND ANTIPARALLEL
		return 2
	return -1

def getBetaAngle ():
	return -1
def getGammaAngle():
	return -1

#NOTE: WE WILL RUN IN BATCHES OF SIMILAR RESOLUTION
pfile = '/Users/fsimon/Google Drive/School/Winter_13/Digital Biology/Project/DigBioProj_One/test.pdb'
appf = pr.parsePDB(pfile, model=1, secondary=True, chain='A', altLoc=False)

# get secondary structure by aParsedPDBfile[i].getSecstr()
nitrox_don = appf.select('element N O') #O can be donor in rare cases, the paper uses this convention
oxygen_acc = appf.select('element O')

for no_d in nitrox_don:
	for ox_a in oxygen_acc:

		#1 Dist(don, acc) < 3.5
		da_dist 	 = pr.calcDistance( no_d , ox_a )
		if( da_dist >= DONOR_ACCEPTOR_MAXDISTANCE ):
			continue

		#2 Dist(h_don, acc) < 3.5
		h_don		 = getHforAtom( appf, no_d ) #Get h_don
		ha_dist 	 = pr.calcDistance(h_don, ox_a)
		if( ha_dist >=  HYDROGEN_ACCEPTOR_MAXDISTANCE):
			continue

		#3 Angle(don, h_don, acc) > 90
		dha_ang = pr.calcAngle( no_d, h_don, ox_a )
		if( dha_ang < DHA_ANGLE ):
			continue

		#4 Angle(don, acc, acc_ante) > 90
		acc_ante 	= getAntecedent(appf, ox_a) #Get acc_ante
		daa_ang 	=  pr.calcAngle( no_d, ox_a, acc_ante)
		if( daa_ang < DAB_ANGLE ):
			continue

		#5 Angle(h_don, acc, acc_ante) > 90
		haa_ang 	= pr.calcAngle( h_don, ox_a, acc_ante )
		if( haa_ang < HAB_ANGLE ):
			continue

		#We have a valid H-Bond with no_d, ox_a, h_don
		print 'HBond elements', h_don.getName(), no_d.getName(), ox_a.getName()
		print 'DA dist', da_dist
		print 'HA dist', ha_dist
		print 'DHA ang', dha_ang
		print 'DAA ang', daa_ang
		print 'HAA ang', haa_ang
		
		beta_ang = getBetaAngle ()
		gamm_ang = getGammaAngle()
		#PUT DATA INTO COLUMN DATA STRUCTURES
		ssindex = getSSIndex(acc_ante)
		if (ssindex == -1):
			#Not a structure we need
			continue
		COLUMN_D_ON  [ssindex].append(da_dist)
		COLUMN_D_OH  [ssindex].append(ha_dist)
		COLUMN_A_NHO [ssindex].append(dha_ang)
		COLUMN_A_HOC [ssindex].append(haa_ang)
		COLUMN_BETA  [ssindex].append(beta_ang)
		COLUMN_GAMMA [ssindex].append(gamm_ang)

COLUMN_D_ON_AV  = [sum(x)/len(x) for x in COLUMN_D_ON]
COLUMN_D_OH_AV  = [sum(x)/len(x) for x in COLUMN_D_OH]
COLUMN_A_NHO_AV = [sum(x)/len(x) for x in COLUMN_A_NHO]
COLUMN_A_HOC_AV = [sum(x)/len(x) for x in COLUMN_A_HOC]
COLUMN_BETA_AV  = [sum(x)/len(x) for x in COLUMN_BETA]
COLUMN_GAMMA_AV = [sum(x)/len(x) for x in COLUMN_GAMMA]

TABLE = [COLUMN_D_ON_AV, COLUMN_D_OH_AV, COLUMN_A_NHO_AV, COLUMN_A_HOC_AV, COLUMN_BETA_AV, COLUMN_GAMMA_AV]
print '        D_ON          D_OH      ANGLE(NHO)    ANGLE(HOC)        BETA         GAMMA   '
print np.array(TABLE).T

'''
for at1 in nitrox:
	for at2 in nitrox:
		if( pr.calcDistance(at1,at2) >= DONOR_ACCEPTOR_MAXDISTANCE ):
			continue
#ASSUME AT1 IS DONOR
		#Assume at1 is Donor and at2 is Acceptor
		h_of_at1 = getHforAtom(appf, at1)
		at1DonorFlag = True
		if( pr.calcDistance(h_of_at1, at2) >=  HYDROGEN_ACCEPTOR_MAXDISTANCE):
			at1DonorFlag = False
			continue
		if ( pr.calcAngle(at1, h_of_at1, at2) <= DHA_ANGLE ):
			at1DonorFlag = False
			continue
		acceptor_antecedent = 
#ASSUME AT2 IS DONOR
		#Assume at2 is Donor and at1 is Acceptor
		h_of_at2 = getHforAtom(appf, at2)
		at2DonorFlag = True
		if( pr.calcDistance(h_of_at2, at1) >  HYDROGEN_ACCEPTOR_MAXDISTANCE):
			at2DonorFlag = False
			continue
'''

'''
def getAntecedent(aParsedPDBfile, anAtom):
	indOfAtom = anAtom.getIndex()
	rng = range(indOfAtom-3, indOfAtom+3) #arbitrary 3 behind and ahead
	rng = [x for x in rng if x>0 & x<aParsedPDBfile.numAtoms()]
	rng.reverse()
	for ind in rng:
		if(aParsedPDBfile[ind].getName() == 'C'):
			return aParsedPDBfile[ind]
'''

'''
# antecedents = 
donors = nitrogens
acceptors = oxygens
for donor in donors:
	for acceptor in acceptors:
		dist = pr.calcDistance(donor, acceptor)
		if(dist >= 3.5):
			continue
		hydrogen = getHforAtom(appf, acceptor)
		hdistance = pr.calcDistance(acceptor, hydrogen)
		if(hdistance >= 2.5):
			continue
		antecedent = getAntecedent(appf, acceptor)
		dha = pr.calcAngle(hydrogen, acceptor, antecedent)
		if(dha <= 90):
			conitnue
		#We have a valid h-bond
		print donor, acceptor
'''