#!/bin/bash

########################## BenchEvalFS2Dpseudo2d #################################
#
# author: Frédéric Darboux <Frederic.Darboux@orleans.inra.fr> (2012-2015)
# version: 1.06.01
# date: 2015-10-29
#
# License Cecill-V2 <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>
#
# (c) CNRS - Universite d'Orleans - INRA (France)
#
# This file is part of FullSWOF_2D software. 
# <https://sourcesup.renater.fr/projects/fullswof-2d/> 
#
# FullSWOF_2D = Full Shallow-Water equations for Overland Flow, 
# in two dimensions of space.
# This software is a computer program whose purpose is to compute
# solutions for 2D Shallow-Water equations.
#
# LICENSE
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# <http://www.cecill.info>.
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##################################################################################

##################################################################################
# The present script "BenchEvalFS2Dpseudo2D" helps in benchmarking the results of
# FullSWOF_2D for pseudo 2D cases. It computes absolute and relative differences
# between a reference result and a result to be evaluated.
##################################################################################

# To do
# - treat cases with different resolutions

VERSION="1.06.01, 2015-10-29"

# Define output flux
STDERR=2

#to have . as decimal separator
LANG=C

if [ $# -ne 5 ]
then
  echo "Usage: BenchEvalFS2Dpseudo2D Reference_filename Ref_HeaderLength Eval_filename Eval_HeaderLength Zmax" >&$STDERR
  echo "Version $VERSION" >&$STDERR
  exit 1;
fi ;

# Define which command for awk
AWK_CMD=gawk
# Define which command for sort
SORT_CMD=sort
# Define which command for cut
CUT_CMD=cut
# Define which command for paste
PASTE_CMD=paste
# Define which command for mktemp
MKTEMP_CMD=mktemp

# Below which absolute value to consider a float as equal to zero
FLOATEQUALZERO=1e-300


#Define which input file to use
REF_FILE=$1
REF_HEADERLENGTH=$2
EVAL_FILE=$3
EVAL_HEADERLENGTH=$4
Z_MAX=$5

#Create temporary files
Ref1Tempfile=$($MKTEMP_CMD Ref1_XXXXXXXX)
Ref2Tempfile=$($MKTEMP_CMD Ref2_XXXXXXXX)
Eval1Tempfile=$($MKTEMP_CMD Eval1_XXXXXXXX)
Eval2Tempfile=$($MKTEMP_CMD Eval2_XXXXXXXX)
Diff1Tempfile=$($MKTEMP_CMD Diff1_XXXXXXXX)
Diff2Tempfile=$($MKTEMP_CMD Diff2_XXXXXXXX)

##function to remove temporary files
function rmtmp(){
rm "$Ref1Tempfile" "$Ref2Tempfile" "$Eval1Tempfile" "$Eval2Tempfile" "$Diff1Tempfile" "$Diff2Tempfile"
return 0
}

##function to compute statistics over diff
function computeStats(){
inputfile=$1
modifiedfile=$($MKTEMP_CMD CMPXXXXXXXX)
#nbdiff=NaN
echo -n "$2" ; $AWK_CMD '$1 == "NaN" {i++} END{printf(" nbdiff==NaN %e\n", i)}' "$inputfile"
#suppress NaN and sort file (sorting will help for min, max and median) for scientific numbers
grep -v 'NaN' "$inputfile" | $SORT_CMD -g > "$modifiedfile"
#nbdiff=0
echo -n "$2" ; $AWK_CMD '$1 == 0 {i++} END{printf(" nbdiff==0 %e\n", i)}' "$modifiedfile"
#nbdiff >0
echo -n "$2" ; $AWK_CMD '$1 > 0 {i++} END{printf(" nbdiff>0 %e\n", i)}' "$modifiedfile"
#nbdiff<0
echo -n "$2" ; $AWK_CMD '$1 < 0 {i++} END{printf(" nbdiff<0 %e\n", i)}' "$modifiedfile"
#min and #max
MINMAX=1
if(test "$(wc -l < "$modifiedfile")" -gt 0); then #at least one line in the file
  MIN=$(head -1 "$modifiedfile")
  echo "$2 min $MIN"
  MAX=$(tail -1 "$modifiedfile")
  echo "$2 max $MAX"
else
  MINMAX=0
  echo "$2 min NaN"
  echo "$2 max NaN"
fi
#mean
echo -n "$2" ; $AWK_CMD '{TOTAL+=$1} END{if(NR>0){printf(" mean %e\n", TOTAL/NR)}else{printf(" mean NaN\n")}}' "$modifiedfile"
#median  (since it is sorted, just take the value at the center)
echo -n "$2" ; $AWK_CMD '{count[NR] = $1} END{if(NR>0){if (NR % 2) {median=count[(NR + 1) / 2];} else {median=(count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0} printf(" median %e\n", median)}else{printf(" median NaN\n")}}' "$modifiedfile"
#Norm L1
echo -n "$2" ; $AWK_CMD '{if ($1 < 0) $1 = -$1; TOTAL+=$1} END{if(NR>0){printf(" L1 %e\n", TOTAL/NR)}else{printf(" L1 NaN\n")}}' "$modifiedfile"
#Norm L2
echo -n "$2" ; $AWK_CMD '{TOTAL+=$1*$1} END{if(NR>0){printf(" L2 %e\n", sqrt(TOTAL/NR))}else{printf(" L2 NaN\n")}}' "$modifiedfile"
#Norm Linfinity
if(test "$MINMAX" -ne 0); then #min and max exists
  echo -n "$2" ; $AWK_CMD 'BEGIN{min='$MIN'; max='$MAX'; if(max<0)max=-max; if(min<0)min=-min; if(min>max)max=min; printf(" Linf %e\n", max)}' "$modifiedfile"
else
  echo "$2 Linf NaN"
fi

rm "$modifiedfile"
return 0
}


##Remove unuseful lines
#Remove header from Reference solution
$AWK_CMD 'FNR>HEADERLENGTH' HEADERLENGTH="$REF_HEADERLENGTH" "$REF_FILE" > "$Ref1Tempfile"
#Remove empty lines
$AWK_CMD '/./' "$Ref1Tempfile" > "$Ref2Tempfile"

#Remove header from Evaluated solution
$AWK_CMD 'FNR>HEADERLENGTH' HEADERLENGTH="$EVAL_HEADERLENGTH" "$EVAL_FILE" > "$Eval1Tempfile"
#Remove empty lines
$AWK_CMD '/./' "$Eval1Tempfile" > "$Eval2Tempfile"

#Compute mean value along each X for Evaluated solution
$AWK_CMD '
{
if(1==FNR){#initialization for first line
prevx=$1
}
x=$1
if (x == prevx){
   if(ZMAX > $7){
      total+=$3
      no++
   }
}else{
   print prevx, "\t", total/no 
   if(ZMAX > $7){
      total=$3
      no=1
   }
   else{
      total=0
      no=0
   }
}
prevx=x
}
END{ #print last computation
   print prevx, "\t", total/no
}
' ZMAX=$Z_MAX "$Eval2Tempfile" > "$Eval1Tempfile"

#check if same number of lines
if (test "$(wc -l < "$Ref1Tempfile")" -ne "$(wc -l < "$Eval1Tempfile")"); then
	echo "Error: Numbers of lines are different." >&"$STDERR";
#	rmtmp;
	exit 1;
fi

#extract relevant columns from Reference File
$CUT_CMD -f1-2 "$Ref2Tempfile"  > "$Ref1Tempfile"

#create a single file
$PASTE_CMD "$Ref1Tempfile" "$Eval1Tempfile" > "$Diff1Tempfile"

#check if Xs are identical
$AWK_CMD '{if (($1 != $3)) {print "Error: X are different at line", FNR, ":", $0; exit 1}}' "$Diff1Tempfile" >&$STDERR
if (test $? -ne 0); then
	echo "Error: Problem while comparing Xs and Ys in the files"  >&$STDERR
	rmtmp;
	exit 1;
fi

#write header in output file

echo    "##############################################################################"
echo    "# Generated by BenchEvalFS2Dpseudo2D version $VERSION"
echo    "# on $(date "+%Y-%m-%d %H:%M:%S") ($(id -un)@$(hostname))"
echo -e "# from reference  file (header = $REF_HEADERLENGTH lines): \t$(pwd)/$REF_FILE"
echo -e "# from evaluation file (header = $EVAL_HEADERLENGTH lines): \t$(pwd)/$EVAL_FILE"
echo    "##############################################################################"

##compute differences and statistics
## for h
#difference in m
$AWK_CMD '{print $4-$2}' "$Diff1Tempfile" > "$Diff2Tempfile"
computeStats "$Diff2Tempfile" DhSI ;
#difference in %
$AWK_CMD '{if($2<-'$FLOATEQUALZERO' || $2>'$FLOATEQUALZERO'){print 100*($4/$2-1)}else{print "NaN"}}' "$Diff1Tempfile" > "$Diff2Tempfile"
computeStats "$Diff2Tempfile" Dh% ;


rmtmp;

exit 0
