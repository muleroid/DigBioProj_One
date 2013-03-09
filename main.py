import prody as pr
import math
import numpy as np
import os
import sheets

DONOR_ACCEPTOR_MAXDISTANCE = 3.5
HYDROGEN_ACCEPTOR_MAXDISTANCE = 2.5
DHA_ANGLE = 90
DAB_ANGLE = 90
HAB_ANGLE = 90

#COLUMNS FOR DATA
#STRUCTURE INDICES ARE: 0=ALPHA, 1=310ALPHA, 2=PARALLEL, 3=ANTIPARALLEL
#COLUMN[STRUCTURE][1-N]
COLUMN_D_ON  = [[],[],[],[]] 
COLUMN_D_OH  = [[],[],[],[]] 
COLUMN_A_NHO = [[],[],[],[]] 
COLUMN_A_HOC = [[],[],[],[]] 
COLUMN_BETA  = [[],[],[],[]] 
COLUMN_GAMMA = [[],[],[],[]] 

def getHforAtom(appf, anAtom):
    resnumneighbors = appf.select('resnum '+str(anAtom.getResnum()))
    curMin = 50
    for at in resnumneighbors:
        if(at.getElement() == 'H'):
            #FILTER OTHER HYDROGENS
            dist = pr.calcDistance(at,anAtom)
            if(dist <= curMin):
            #print pr.calcDistance(at,anAtom), at.getResnum(), at.getName()
                curMin = dist
                target = at
    return at

def setHcoords(nAtom, oAtom, cAtom):
    # get coordinates
    nCoords = nAtom.getCoords()
    oCoords = oAtom.getCoords()
    cCoords = cAtom.getCoords()
    # C-O vector
    co = np.subtract(cCoords,oCoords)
    adjustment = np.divide(co,np.linalg.norm(co))
    hCoords = np.add(nCoords,adjustment)
    return hCoords

def getAntecedent (appf, anAtom):
    aminoGroup = appf.select('resnum ' + str(anAtom.getResnum()))
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

def getBetaAngle(appf,cAtom,oAtom,hAtom):
    # first determine the nAtom
    aminoGroup = appf.select('resnum ' + str(cAtom.getResnum()))
    for at in aminoGroup:
        if(at.getName() == 'N'):
            nAtom = at
            # get coordinates
            cCoords = cAtom.getCoords()
            oCoords = oAtom.getCoords()
            hCoords = hAtom.getCoords()
            nCoords = nAtom.getCoords()
        # get relevant vectors
        oc = np.subtract(oCoords,cCoords)
        nc = np.subtract(nCoords,cCoords)
        ho = np.subtract(hCoords,oCoords)
        # find norm vector to N-C-O plane
        n1 = np.cross(oc,nc)
        n1_unit = np.divide(n1,np.linalg.norm(n1))
        # find projection of H-O onto plane
        out = np.dot(ho,n1_unit)
        ho_ip = np.subtract(ho,np.multiply(n1_unit,out))
        ang = np.linalg.norm(ho_ip)/np.linalg.norm(ho)
        ang = math.acos(ang)
        # we can think of H-O vector as the normal vector to a plane
        #ang = math.acos(np.dot(ho,n1)/(np.linalg.norm(ho)*np.linalg.norm(n1)))
        ang = ang*180/math.pi
        # then we just need to adjust output by 90 degrees to get correct answer
        return ang

def getGammaAngle(appf,cAtom,oAtom,hAtom):
    # first determine the nAtom
    aminoGroup = appf.select('resnum ' + str(cAtom.getResnum()))
    for at in aminoGroup:
        if(at.getName() == 'N'):
            nAtom = at
        # get coordinates
    cCoords = cAtom.getCoords()
    oCoords = oAtom.getCoords()
    hCoords = hAtom.getCoords()
    nCoords = nAtom.getCoords()
    # get necessary vectors
    oc = np.subtract(oCoords,cCoords)
    nc = np.subtract(nCoords,cCoords)
    ho = np.subtract(hCoords,oCoords)
    n1 = np.cross(oc,nc)
    n1_unit = np.divide(n1,np.linalg.norm(n1))
    # get projection of H-O in O-C direction
    oc_unit = np.divide(oc,np.linalg.norm(oc))
    #print oc_unit
    hproj = np.dot(ho,oc_unit)
    # get projection of H-O onto N-C-O plane
    out = np.dot(ho,n1_unit)
    ho_ip = np.subtract(ho,np.multiply(n1_unit,out))
    ang = hproj/np.linalg.norm(ho_ip)
    ang = math.acos(ang)
    ang = ang*180/math.pi
    return ang*-1
    
# find the number of models in a file
def getNumMdls(pfile):
    fp = open(pfile,'r')
    models = 0
    for line in fp:
        #print line[0:4]
        if(line[0:5] == "MODEL"):
            models = models+1
    fp.close()
    if(models == 0):
        return 1
    return models

#NOTE: WE WILL RUN IN BATCHES OF SIMILAR RESOLUTION
#pfile = '/Users/fsimon/Google Drive/School/Winter_13/Digital Biology/Project/DigBioProj_One/test.pdb'
#pfile = '1A6S_A_H.pdb'

def checkHBonds(appf,no_d, ox_a, h_don, acc_ante, ssindex):
    #1 Dist(don, acc) < 3.5
    da_dist = pr.calcDistance( no_d , ox_a )
    if( da_dist >= DONOR_ACCEPTOR_MAXDISTANCE ):
        return

    #2 Dist(h_don, acc) < 3.5
    ha_dist  = pr.calcDistance(h_don, ox_a)
    if( ha_dist >=  HYDROGEN_ACCEPTOR_MAXDISTANCE):
        return
            
    #3 Angle(don, h_don, acc) > 90
    dha_ang = pr.calcAngle( no_d, h_don, ox_a )
    if( dha_ang < DHA_ANGLE ):
        return
            
    #4 Angle(don, acc, acc_ante) > 90
    #daa_ang =  pr.calcAngle( no_d, ox_a, acc_ante)
    #if( daa_ang < DAB_ANGLE ):
    #    return
            
    #5 Angle(h_don, acc, acc_ante) > 90
    haa_ang = pr.calcAngle( h_don, ox_a, acc_ante )
    if( haa_ang < HAB_ANGLE ):
        return

    #We have a valid H-Bond with no_d, ox_a, h_don
    #print 'HBond elements', h_don.getName(), no_d.getName(), ox_a.getName()
    #print 'DA dist', da_dist
    #print 'HA dist', ha_dist
    #print 'DHA ang', dha_ang
    #print 'DAA ang', daa_ang
    #print 'HAA ang', haa_ang
            
    beta_ang = getBetaAngle (appf,acc_ante,ox_a,h_don)
    gamm_ang = getGammaAngle(appf,acc_ante,ox_a,h_don)
        
    #PUT DATA INTO COLUMN DATA STRUCTURES
#    print ssindex
    #towrite = ",".join([str(h_don.getResnum()),str(h_don.getCoords())])
    #fp.write(towrite+"\r\n")
    COLUMN_D_ON  [ssindex].append(da_dist)
    COLUMN_D_OH  [ssindex].append(ha_dist)
    COLUMN_A_NHO [ssindex].append(dha_ang)
    COLUMN_A_HOC [ssindex].append(haa_ang)
    COLUMN_BETA  [ssindex].append(beta_ang)
    COLUMN_GAMMA [ssindex].append(gamm_ang)
    return

def parseHelices(appf):
    # get secondary structure by aParsedPDBfile[i].getSecstr()
    nitrox_don = appf.select('element N O') #O can be donor in rare cases, the paper uses this convention
    oxygen_acc = appf.select('element O')
    for no_d in nitrox_don:
        n_secStruct = getSSIndex(no_d)
        # we're only concerned with helices
        if(n_secStruct < 0 or n_secStruct == 2):
            continue
        for ox_a in oxygen_acc:
            o_secStruct = getSSIndex(ox_a)
            if(o_secStruct != n_secStruct):
                continue
            acc_ante = getAntecedent(appf, ox_a) #Get acc_ante
            ssindex = getSSIndex(acc_ante) # get ssindex
            h_don = getHforAtom(appf,no_d)
            #newCoords = setHcoords(no_d,ox_a,acc_ante)
            #h_don.setCoords(newCoords)
            checkHBonds(appf,no_d,ox_a,h_don,acc_ante,ssindex)
    return

# do the sheets
def parseSheets(appf, los):
    strands = sheets.buildStrands(appf, los)
    numStrands = len(strands)
    for i in np.arange(numStrands-2):
        base = strands[i]
        comp = strands[i+1]
        if(comp.sense == 0):
            continue
        for no_d in base.getNAtoms():
            for ox_a in comp.getOAtoms():
                acc_ante = getAntecedent(appf, ox_a)
                if(comp.getSense() > 0):
                    ssindex = 2
                else:
                    ssindex = 3
                h_don = getHforAtom(appf,no_d)
                newCoords = setHcoords(no_d,ox_a,acc_ante)
                h_don.setCoords(newCoords)
                checkHBonds(appf,no_d,ox_a,h_don,acc_ante,ssindex)
    return

# main function
def runThrough(pfile):
    # initial setup
    print "Running through " + pfile + "..."
    numMdls = getNumMdls(pfile)
    #print numMdls
    appf = pr.parsePDB(pfile, model=numMdls, secondary=True, chain='A', altLoc=False)
    los = sheets.initializeList(pfile)
    parseHelices(appf)
    parseSheets(appf,los)

    # print out our results
    COLUMN_D_ON_AV  = [sum(x)/len(x) for x in COLUMN_D_ON if len(x) > 0]
    COLUMN_D_OH_AV  = [sum(x)/len(x) for x in COLUMN_D_OH if len(x) > 0]
    COLUMN_A_NHO_AV = [sum(x)/len(x) for x in COLUMN_A_NHO if len(x) > 0]
    COLUMN_A_HOC_AV = [sum(x)/len(x) for x in COLUMN_A_HOC if len(x) > 0]
    COLUMN_BETA_AV  = [sum(x)/len(x) for x in COLUMN_BETA if len(x) > 0]
    COLUMN_GAMMA_AV = [sum(x)/len(x) for x in COLUMN_GAMMA if len(x) > 0]

    TABLE = [COLUMN_D_ON_AV, COLUMN_D_OH_AV, COLUMN_A_NHO_AV, COLUMN_A_HOC_AV, COLUMN_BETA_AV, COLUMN_GAMMA_AV]
    print '        D_ON          D_OH      ANGLE(NHO)    ANGLE(HOC)        BETA         GAMMA   '
    print np.array(TABLE).T

# code to run through all the files and stuff
AGGR_D_ON  = [0,0,0,0] 
AGGR_D_OH  = [0,0,0,0] 
AGGR_A_NHO = [0,0,0,0] 
AGGR_A_HOC = [0,0,0,0] 
AGGR_BETA  = [0,0,0,0] 
AGGR_GAMMA = [0,0,0,0]

LEN_D_ON = [0,0,0,0]
LEN_D_OH = [0,0,0,0]
LEN_A_NHO = [0,0,0,0]
LEN_A_HOC = [0,0,0,0]
LEN_BETA = [0,0,0,0]
LEN_GAMMA = [0,0,0,0]

pFiles = os.listdir(os.getcwd())
for pfile in pFiles:
    # reset everything
    COLUMN_D_ON  = [[],[],[],[]] 
    COLUMN_D_OH  = [[],[],[],[]] 
    COLUMN_A_NHO = [[],[],[],[]] 
    COLUMN_A_HOC = [[],[],[],[]] 
    COLUMN_BETA  = [[],[],[],[]] 
    COLUMN_GAMMA = [[],[],[],[]] 
    runThrough(pfile)
    # aggregate that shit
    AGGR_D_ON = [(AGGR_D_ON[i] + sum(COLUMN_D_ON[i])) for i in range(4)]
    AGGR_D_OH = [(AGGR_D_OH[i] + sum(COLUMN_D_OH[i])) for i in range(4)]
    AGGR_A_NHO = [(AGGR_A_NHO[i] + sum(COLUMN_A_NHO[i])) for i in range(4)]
    AGGR_A_HOC = [(AGGR_A_HOC[i] + sum(COLUMN_A_HOC[i])) for i in range(4)]
    AGGR_BETA = [(AGGR_BETA[i] + sum(COLUMN_BETA[i])) for i in range(4)]
    AGGR_GAMMA = [(AGGR_GAMMA[i] + sum(COLUMN_GAMMA[i])) for i in range(4)]

    LEN_D_ON = [(LEN_D_ON[i] + len(COLUMN_D_ON[i])) for i in range(4)]
    LEN_D_OH = [(LEN_D_OH[i] + len(COLUMN_D_OH[i])) for i in range(4)]
    LEN_A_NHO = [(LEN_A_NHO[i] + len(COLUMN_A_NHO[i])) for i in range(4)]
    LEN_A_HOC = [(LEN_A_HOC[i] + len(COLUMN_A_HOC[i])) for i in range(4)]
    LEN_BETA = [(LEN_BETA + len(COLUMN_BETA[i])) for i in range(4)]
    LEN_GAMMA = [(LEN_GAMMA + len(COLUMN_GAMMA[i])) for i in range(4)]

# calculate the averages
AVG_D_ON = [(AGGR_D_ON[i]/LEN_D_ON[i]) for i in range(4)]
AVG_D_OH = [(AGGR_D_OH[i]/LEN_D_OH[i]) for i in range(4)]
AVG_A_NHO = [(AGGR_A_NHO[i]/LEN_A_NHO[i]) for i in range(4)]
AVG_A_HOC = [(AGGR_A_HOC[i]/LEN_A_HOC[i]) for i in range(4)]
AVG_BETA = [(AGGR_BETA[i]/LEN_BETA[i]) for i in range(4)]
AVG_GAMMA = [(AGGR_GAMMA[i]/LEN_GAMMA[i]) for i in range(4)]
TABLE = [AVG_D_ON,AVG_D_OH,AVG_A_NHO,AVG_A_HOC,AVG_BETA,AVG_GAMMA]
print np.array(TABLE).T
