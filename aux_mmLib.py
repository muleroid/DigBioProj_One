#############################################
# Contains auxiliary functions and code to
# run the examples from mmLib
#
#############################################

############
# To run examples
# import subprocess
# paex = '/Users/fsimon/Downloads/pymmlib-1.2.0/examples/' #shortcut path
# subprocess.call([paex+'pdb_parser.py', '/Users/fsimon/Desktop/nrpdb/1A2Z_A.pdb'])
# 
# NOTES
# 1. ITERATE OVER SECONDARY STRUCTURES
#	 struct = FileIO.LoadStructure(file = path)
#    for ahelix in struct.iter_alpha_helicies():
#	 for bsheet in struct.iter_beta_sheets():
#	 for site in struct.iter_sites():
# 
# 2. ITERATING IN BATCHES
#
#	def neighborhood(iterable):
#	    iterator = iter(iterable)
#	    prev = None
#	    item = iterator.next()  # throws StopIteration if empty.
#	    for next in iterator:
#	        yield (prev,item,next)
#	        prev = item
#	        item = next
#	    yield (prev,item,None)
#	
#	for prev,item,next in neighborhood(l):
#	    print prev, item, next
# 
# 3. AUXILIARY FUNCTIONS
# 	calc_distance(a1, a2)
#   calc_torsion_angle(a1, a2, a3, a4)
#        Calculates the torsion angle between the four argument atoms.
#   
#
############

import os, sys
import numpy as np
from mmLib import FileIO
from mmLib import AtomMath
from itertools import islice, chain


#path = '/Users/fsimon/Desktop/nrpdb/1A2Z_A.pdb'
#protein = open(path, "r")

#struct = FileIO.LoadStructure(file = path)

def neighborhood(iterable):
    iterator = iter(iterable)
    prev = None
    item = iterator.next()  # throws StopIteration if empty.
    for next in iterator:
        yield (prev,item,next)
        prev = item
        item = next
    yield (prev,item,None)

for prev,item,next in neighborhood(struct.iter_all_atoms()):
	rad = AtomMath.calc_angle(prev, item, next)
	print prev, item, next
	if(rad != None):
		print np.degrees(rad)

