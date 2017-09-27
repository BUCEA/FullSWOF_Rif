#!/bin/bash

########################## reportDiffFS2D #########################################
#
# author: Frédéric Darboux <Frederic.Darboux@orleans.inra.fr> (2012-2015)
# version: 1.06.00
# date: 2015-02-18
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
# The present script "reportDiffFS2D" helps in benchmarking the results of
# FullSWOF_2D. It displays the biggest absolute and relative differences in the
# file created by the script difCMPFS2D.sh
##################################################################################

VERSION="1.06.00, 2015-02-18"

# Define output flux
STDERR=2

#to have . as decimal separator
LANG=C

if [ $# -ne 1 ]
then
  echo "Usage: reportDiffFS2D diff_FS2D_filename" >&$STDERR
  echo "Version $VERSION" >&$STDERR
  exit 1;
fi ;

# Define which command for awk
AWK_CMD=gawk

# Below which absolute value to consider a float as equal to zero
FLOATEQUALZERO=1e-300

# number of header lines in the input files
HEADERLENGTH=7

#Define which input file to use
DIFF_FILE=$1
 
$AWK_CMD '
NR<=HEADERLENGTH {next}
BEGIN{
minabs=0;
maxabs=0;
minrel=0;
maxrel=0;
}
{
if($3!="NaN" && $3!=0){
	if($3<minabs) minabs=$3;
	if($3>maxabs) maxabs=$3;
}
if($4!="NaN" && $4!=0){
	if($4<minrel) minrel=$4;
	if($4>maxrel) maxrel=$4;
}
}
END{
if(minabs<0)minabs=-minabs;
if(minabs>maxabs)maxabs=minabs;
if(minrel<0)minrel=-minrel;
if(minrel>maxrel)maxrel=minrel;
if(maxabs<FLOATEQUALZERO && maxrel<FLOATEQUALZERO) 
	printf("Results are identical.\n")
else
	printf("Maximum differences are: %e (absolute) and %e (relative).\n", maxabs, maxrel)
}
' HEADERLENGTH=$HEADERLENGTH FLOATEQUALZERO=$FLOATEQUALZERO  "$DIFF_FILE"

exit 0
