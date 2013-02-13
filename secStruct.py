#! /usr/bin/env python

# parses information about secondary structures
import os
import math

def parseHelix(pdbFile):
    helixBounds = []
    helix310Bounds = []
    for pdbRecord in pdbFile:
        if(len(pdbRecord) >= 8):
            if(pdbRecord[0:5] == "SHEET"):
                beg = int(pdbRecord[22:26].strip())
                end = int(pdbRecord[33:37].strip())
                if(int(pdbRecord[38:40]) == 5):
                    helix310Bounds.append(beg,end)
                else:
                    helixBounds.append((beg,end))
    return helixBounds
    
def parseSheet(pdbFile):
    # need to ask prof Scott about anti-parallel vs. parallel
    parSheetBounds = []
    antiparSheetBounds = []
    return 0
