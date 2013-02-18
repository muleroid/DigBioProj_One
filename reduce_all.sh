#! /bin/bash
for i in $(ls *.pdb); do
    b=$(basename $i .pdb)
    output="${b}_H.pdb"
    reduce -build $i > $output
done