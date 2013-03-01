import prody as pr
import numpy as np

class Strand():
    sense = 0
    donors = []
    acceptors = []
    def __init__(self, sense, dons, accs):
        self.sense = sense
        self.donors = dons
        self.acceptors = accs

    def getNAtoms(self):
        return self.donors

    def getOAtoms(self):
        return self.acceptors

    def getSense(self):
        return self.sense

#Get all beta sheets, separate between parallel and antiparallel
#Get The Whole sheet

def initializeList(pfile):
    #filepath = "/Users/fsimon/Google Drive/School/Winter_13/Digital Biology/Project/DigBioProj_One/test.pdb"
    fp = open(pfile,'r')
    list_of_strands = []

    strand_unique_id = 0
    for line in fp:
        if("SHEET" in line[0:5]):       
            strand_number     = int(line[7:10])
            sheet_identifier  = line[11:14]
            number_of_strands = int(line[14:16]) #numstrands in current sheet
            chain_identifier  = line[21]
            start_res         = int(line[22:26])
            chain_identifier2 = line[32]
            end_res           = int(line[33:37])
            strand_sense      = int(line[38:40])
            if(strand_sense == 0):
                #We're on a new strand
                strand_unique_id  = strand_unique_id + 1
            list_of_strands.append([strand_sense, start_res, end_res, strand_unique_id])
        if("LINK" in line[0:4]):
            break

    fp.close()
    print list_of_strands
    return list_of_strands

# pull necessary atoms for each strand
def buildStrands(appf, list_of_strands):
    curDons = []
    curAccs = []
    strands = []
    for s in list_of_strands:
        sense = s[0]
        start = s[1]
        end = s[2]
        for resNum in np.arange(start,end+1):
            curGroup = appf.select('resnum ' + str(resNum))
            for atom in curGroup:
                if(atom.getName() == 'N'):
                    curDons.append(atom)
                if(atom.getName() == 'O'):
                    curAccs.append(atom)
        strands.append(Strand(sense, curDons, curAccs))
        curDons = []
        curAccs = []
    return strands

def test(pfile):
    appf = pr.parsePDB(pfile, model=1, secondary=True, chain='A', altLoc=False)
    los = initializeList(pfile)
    strands = buildStrands(appf, los)
    # print out strands
    for strand in strands:
        dons = strand.getNAtoms()
        accs = strand.getOAtoms()
        for n in dons:
            print n.getResnum()
            print n.getResname()
    return

test('1A2Z_A.pdb')

##def getUniqueIDofStrand (residueNumber, listOfStrands):
##    for strand in listOfStrands:
##        if( (residueNumber >= strand[1]) & (residueNumber <= strand[2]) ):
##            return strand[3]
##
##def returnParOrAnti(residueNumber, listOfStrands):
##    thisUID = getUniqueIDofStrand(residueNumber, listOfStrands)
##    parallel = None
##
##    for i in range(0, len(listOfStrands)):
##        #Search list to find correct strand.
##        strand = listOfStrands[i]
##        if( (residueNumber >= strand[1]) & (residueNumber <= strand[2]) ):
##            #Found the strand we want.
##
##            if(strand[0] == 0):
##                return "BASE"
##
##            for strand in listOfStrands:
##                #Toggle parallelbool from origin to sink to figure out anti vs parallel
##                if(thisUID == strand[3]):
##                    if(residueNumber > strand[2]): #went too far
##                        break
##                    #Only toggle in the same sheet.
##                    if(strand[0] == -1):
##                        print strand[1], strand[2]
##                        if(parallel == None):
##                            parallel = False
##                            continue
##                        parallel = not parallel
##                    if(strand[0] == 1):
##                        if(parallel == None):
##                            parallel = True
##                            continue
##                        parallel = parallel
##
##            if(parallel == True):
##                return "PARALLEL"
##            else:
##                return "ANTIPARALLEL"
##
##
##                      
##                #look back until we found a 0 sense strand.
##
##thelist = initializeList()
##print returnParOrAnti(71, thelist)
##
##
##                #initial_resname   = line[17:20]
##                        #code_for_insert   = line[27]
##        #termin_resid_name = line[28:31]
##                #code_for_insert2  = line[37]
