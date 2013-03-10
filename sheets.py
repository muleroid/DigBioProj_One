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
#    print list_of_strands
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
            if(curGroup == None):
                continue
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
            print n.getSecstr()
    return
