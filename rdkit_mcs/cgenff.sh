#!/bin/bash

export cgenff=/Data/cbu/cbureuterfs/software/silcsbio.2026.1-alpha/cgenff/cgenff

while read i
do

#$cgenff -v -a $i.mol2 > $i.str
#sed -n -e '/Toppar/,/END/ p' $i.str > $i.rtf
#sed -n -e '/flex/,/RETURN/ p' $i.str > $i.prm

./regroup.awk $i.str > "$i"_rg.str
./regroup.awk $i.rtf > "$i"_rg.rtf

done < <(cat mol_list.txt)
