#! /usr/env/bin python

import os, re, shutil, StringIO
nrpdbpath = "./"
directory = nrpdbpath
dirbuck0  = nrpdbpath + "buck0/"
dirbuck1  = nrpdbpath + "buck1/"
dirbuck2  = nrpdbpath + "buck2/"
dirbuck3  = nrpdbpath + "buck3/"
dirbuck4  = nrpdbpath + "buck4/"
dirbuck5  = nrpdbpath + "buck5/"
dirbuck6  = nrpdbpath + "buck6/"
dirbuck7  = nrpdbpath + "buck7/"
dirbuck8  = nrpdbpath + "buck8/"
dirbuck_nores  = nrpdbpath + "buck_nores/"


allresolutions = []

int1 = 0
for item in os.listdir(directory):
    int1 = int1+1
    print int1
    resolution = -2
    if(item.startswith(".")):
        continue
    if(not item.endswith("pdb")):
        continue
    fil = open(directory+item)
    file_lines = fil.readlines()
    for line in file_lines:
        if ("RESOLUTION" in line):
            if ("2 RESOLUTION. NOT APPLICABLE." in line):
                resolution = -1
                break
            m = re.search("([0-9]*.[0-9]+)\s+ANGSTROMS.", line)
            if m:
                resolution = float(m.groups()[0])
				#print resolution
                break

    fil.close()
    allresolutions.append(float(resolution))
	#Now we have the resolution.
    if(resolution == -1):
        print "YO"
        os.system("mv " + item + " " + dirbuck_nores)
            #wrfi = open(dirbuck_nores+item, 'w')
            #wrfi.writelines(file_lines)
            #wrfi.close()
            #fil.close()
    elif((resolution > 0.5)) & ((resolution <= 1.0)):
        os.system("mv " + item + " " + dirbuck0)
#wrfi = open(dirbuck0+item, 'w')
            #wrfi.writelines(file_lines)
            #wrfi.close()
            #fil.close()
    elif((resolution > 1)) & ((resolution <= 1.25)):
        os.system("mv " + item + " " + dirbuck1)
#wrfi = open(dirbuck1+item, 'w')
            #wrfi.writelines(file_lines)
            #wrfi.close()
            #fil.close()	
    elif((resolution > 1.25) & (resolution <= 1.5)):
        os.system("mv " + item + " " + dirbuck2)
#wrfi = open(dirbuck2+item, 'w')
            #wrfi.writelines(file_lines)
            #wrfi.close()
            #fil.close()	
    elif((resolution > 1.5) & (resolution <= 1.75)):
        os.system("mv " + item + " " + dirbuck3)
            #wrfi = open(dirbuck3+item, 'w')
            #wrfi.writelines(file_lines)
            #wrfi.close()
            #fil.close()	
    elif((resolution > 1.75) & (resolution <= 2)):
        os.system("mv " + item + " " + dirbuck4)
            #wrfi = open(dirbuck4+item, 'w')
            #wrfi.writelines(file_lines)
            #wrfi.close()
            #fil.close()	
    elif((resolution > 2) & (resolution <= 2.25)):
        os.system("mv " + item + " " + dirbuck5)
            #wrfi = open(dirbuck5+item, 'w')
            #wrfi.writelines(file_lines)
            #wrfi.close()
            #fil.close()
    elif((resolution > 2.25) & (resolution <= 2.5)):
        os.system("mv " + item + " " + dirbuck6)
#wrfi = open(dirbuck6+item, 'w')
            #wrfi.writelines(file_lines)
            #wrfi.close()
            #fil.close()	
    elif((resolution > 2.5) & (resolution <= 2.75)):
        os.system("mv " + item + " " + dirbuck7)
            #wrfi = open(dirbuck7+item, 'w')
            #wrfi.writelines(file_lines)
            #wrfi.close()
            #fil.close()	
    elif((resolution > 2.75)):
        os.system("mv " + item + " " + dirbuck8)
            #wrfi = open(dirbuck8+item, 'w')
            #wrfi.writelines(file_lines)
            #wrfi.close()
            #fil.close()
    else:
        print "ERROR in " + item
