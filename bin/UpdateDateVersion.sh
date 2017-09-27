#!/bin/bash

########################## UpdateDateVersion #########################################
#
# author: Frédéric Darboux <Frederic.Darboux@orleans.inra.fr> (2012-2015)
# version: 1.06.00
# date: 2015-04-03
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
# The present script "UpdateDateVersion" helps in maintaining FullSWOF_2D
# by updating consistently the version number and date of FullSWOF_2D
# in all relevant files, and compiling the documentation files
##################################################################################

SCRIPTVERSION="1.07.00, 2015-04-03."

# Define output flux
STDERR=2

#to have . as decimal separator
LANG=C

if [ $# -ne 2 ]
then
  echo "Usage: UpdateDateVersion YYYY-MM-DD X.YY.ZZ" >&$STDERR
  echo "Version $SCRIPTVERSION" >&$STDERR
  exit 1;
fi ;

# Define which command for sed
SED_CMD=sed
# Define which command for mktemp
MKTEMP_CMD=mktemp
# Define which command for html2pdf
EGREP_CMD=egrep
# Define which command for doxygen
DOXYGEN_CMD=doxygen

#Read input values
DATE=$1
VERSION=$2

### Check for programs
ALLHERE=0
command -v $SED_CMD >/dev/null 2>&1 || { echo >&2 "Error: program" $SED_CMD "required but not installed"; ALLHERE=1; }
command -v $MKTEMP_CMD >/dev/null 2>&1 || { echo >&2 "Error: program" $MKTEMP_CMD "required but not installed"; ALLHERE=1; }
command -v $EGREP_CMD >/dev/null 2>&1 || { echo >&2 "Error: program" $EGREP_CMD "required but not installed"; ALLHERE=1; }
command -v $DOXYGEN_CMD >/dev/null 2>&1 || { echo >&2 "Error: program" $DOXYGEN_CMD "required but not installed"; ALLHERE=1; }
if [ $ALLHERE -eq 1 ] #at least one program missing...
then
 exit 1 ;
fi

###Update README.txt
if [ -f README.txt ] #does the file exists?
then
  #it exists ==> update
  STRING='#### version '$VERSION', '$DATE'' ;
  Tempfile=$($MKTEMP_CMD Tempfile_XXXXXXXX) ;
  $SED_CMD "2s/.*/$STRING/" README.txt > "$Tempfile" && mv -f "$Tempfile" README.txt;
else
  #it does not exist ==> error
  echo "Error: file README.txt is not in the current directory" >&$STDERR ;
  exit 1;
fi

###Update misc.hpp
if [ -f ./Headers/libparameters/misc.hpp ] #does the file exists?
then
  #it exists ==> update
  STRING="#define VERSION \"FullSWOF_2D version $VERSION, $DATE\"";
  Tempfile=$($MKTEMP_CMD Tempfile_XXXXXXXX);
  $SED_CMD "1,/^#define VERSION.*/s/^#define VERSION.*/$STRING/" ./Headers/libparameters/misc.hpp > "$Tempfile" && mv -f "$Tempfile" ./Headers/libparameters/misc.hpp;

  STRING=" * @version $VERSION";
  Tempfile=$($MKTEMP_CMD Tempfile_XXXXXXXX);
  $SED_CMD "1,/^ \* @version .*/s/^ \* @version .*/$STRING/" ./Headers/libparameters/misc.hpp > "$Tempfile" && mv -f "$Tempfile" ./Headers/libparameters/misc.hpp;

  STRING=" * @date $DATE";
  Tempfile=$($MKTEMP_CMD Tempfile_XXXXXXXX);
  $SED_CMD "1,/^ \* @date .*/s/^ \* @date .*/$STRING/" ./Headers/libparameters/misc.hpp > "$Tempfile" && mv -f "$Tempfile" ./Headers/libparameters/misc.hpp;
else
  #it does not exist ==> error
  echo "Error: file misc.hpp is not in the ./Headers/libparameters/ directory" >&$STDERR ;
  exit 1;
fi

###Update Documentation.tex
if [ -f ./doc/Documentation.tex ] #does the file exists?
then
  #it exists ==> update
  STRING="\\\newcommand{\\\version}{$VERSION ($DATE)}";
  Tempfile=$($MKTEMP_CMD Tempfile_XXXXXXXX);
  $SED_CMD "1,/^\\\newcommand{\\\version}{.*/s/^\\\newcommand{\\\version}{.*/$STRING/" ./doc/Documentation.tex > "$Tempfile" && mv -f "$Tempfile" ./doc/Documentation.tex;

  #Compilation script from Doxygen makefile
  cd doc|| exit 1
  rm -f Documentation.aux Documentation.log Documentation.blg Documentation.bbl Documentation.out Documentation.toc Documentation.pdf
  pdflatex Documentation
  bibtex Documentation
  pdflatex Documentation
  latex_count=5 ; 
  while $EGREP_CMD -s 'Rerun (LaTeX|to get cross-references right)' Documentation.log && [ $latex_count -gt 0 ] ;
    do 
      echo "Rerunning latex...." ;
      pdflatex Documentation ;
      latex_count=$(expr $latex_count - 1);
    done
  rm -f Documentation.aux Documentation.log Documentation.blg Documentation.bbl Documentation.out Documentation.toc
  cd .. || exit 1
else
  #it does not exist ==> error
  echo "Error: file Documentation.tex is not in the ./doc directory" >&$STDERR ;
  exit 1;
fi

###Update Doxygen_config_file_latex
FILEFOUND=yes
if [ -f ./doc/Doxygen_config_file_latex ] #does the file exists?
then
  #it exists ==> update
  STRING="PROJECT_NUMBER         = \"v$VERSION\ ($DATE)\"";
  Tempfile=$($MKTEMP_CMD Tempfile_XXXXXXXX);
  $SED_CMD "1,/^PROJECT_NUMBER         =.*/s/^PROJECT_NUMBER         =.*/$STRING/" ./doc/Doxygen_config_file_latex > "$Tempfile" && mv -f "$Tempfile" ./doc/Doxygen_config_file_latex;
else
  #it does not exist ==> error
  echo "Error: file Doxygen_config_file_latex is not in the ./doc directory" >&$STDERR ;
  FILEFOUND=no
  exit 1;
fi

###Update Doxygen_config_file_html
if [ -f ./doc/Doxygen_config_file_html ] #does the file exists?
then
  #it exists ==> update
  STRING="PROJECT_NUMBER         = \"v$VERSION ($DATE)\"";
  Tempfile=$($MKTEMP_CMD Tempfile_XXXXXXXX);
  $SED_CMD "1,/^PROJECT_NUMBER         =.*/s/^PROJECT_NUMBER         =.*/$STRING/" ./doc/Doxygen_config_file_html > "$Tempfile" && mv -f "$Tempfile" ./doc/Doxygen_config_file_html;
else
  #it does not exist ==> error
  echo "Error: file Doxygen_config_file_html is not in the ./doc directory" >&$STDERR ;
  FILEFOUND=no
  exit 1;
fi

###Compile Doxygen
if [[ "$FILEFOUND" == yes ]] #doxygen config files found?
then
  #compilation doxygen LaTeX
  rm -rf doc/latex/
  $DOXYGEN_CMD doc/Doxygen_config_file_latex
  cd doc/latex/ || exit 1
  make
  cp refman.pdf ../
  cd ../../ || exit 1
  rm -rf doc/latex/

  #compilation doxygen html
  rm -rf doc/html/ 
  $DOXYGEN_CMD doc/Doxygen_config_file_html
fi

exit 0
