#!/bin/bash

## LigParGen run on one of the PC-lab systems
## pga043@129.177.190.61 
## conda activate boss

export BOSSdir=/Home/siv32/pga043/Software/boss

while read i
do

echo "$i"

## cgen = CM1A-LBCC working with BOSS5.0

ligpargen -i ../"$i".mol2 -cgen CM1A-LBCC -c 1 -cgenb CM1A-LBCC -cb 1

rm *.q.* *gmx* *tinker* *openmm* *desmond* *xplor* *lammps*

#python replaceAtomTypes.py "$i"

#cp "$i".tmp.rtf ../"$i".rtf
#cp "$i".tmp.prm ../"$i".prm

python opls2charmm.py "$i"

done < <(cat ../mol_list.txt) 

## ligpargen: http://zarbi.chem.yale.edu/ligpargen/index.html

## charge model:  1.14*CM1A-LBCC

##Notes:
##   CM1A is automatically scaled by 1.14 in neutral molecules.
##   CM1A-LBCC is just for neutral molecules and it is also scaled by 1.14.

