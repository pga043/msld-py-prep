#! /bin/bash

export MMTSBDIR=/Home/siv32/pga043/Software/toolset
export PATH=$PATH:$MMTSBDIR/perl:$MMTSBDIR/bin

cutoff=$1
if [[ $cutoff == '' ]] ; then
   cutoff=12
fi
echo "Using a solvate cutoff of $cutoff"

name=`more build*/name`

InDir=./build.${name}
OutDir=solv_prep
Dir=`pwd`
LgLig=large_lig.pdb

cp -r $InDir ./$OutDir
convpdb.pl -solvate -cubic -cutoff $cutoff -ions CLA:1 -removeclashes $InDir/$LgLig | grep TIP3 | convpdb.pl -segnames > $OutDir/solvent.pdb
convpdb.pl -solvate -cubic -cutoff $cutoff -ions CLA:1 -removeclashes $InDir/$LgLig | grep LG | convpdb.pl -segnames | sed "s/PRO0/LIG/g" > $OutDir/LgLig.pdb

convpdb.pl -solvate -cubic -cutoff $cutoff -ions CLA:1 -removeclashes $InDir/$LgLig | grep CLA | convpdb.pl -segnames | sed "s/HETATM/ATOM  /g" > $OutDir/ions.pdb

XYZ=`awk 'BEGIN {CC=0} {if ($1=="ATOM" && $2==1) {if (CC==0) {CC=1; X=$6; Y=$7; Z=$8} else {print X-$6,Y-$7,Z-$8}}}' $OutDir/LgLig.pdb $InDir/$LgLig`

echo $XYZ


cd $OutDir
FILES=`ls 1.pdb 2.pdb`
cd $Dir

for FILE in $FILES
do

  convpdb.pl -translate $XYZ -segnames $InDir/$FILE | sed "s/PRO0/LIG/g" > $OutDir/$FILE

done
