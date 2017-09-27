########################## Makefile ##############################################
#
# author: Ulrich Razafison <ulrich.razafison@math.cnrs.fr> (2008)
# author: Christian Laguerre <christian.laguerre@math.cnrs.fr> (2010)
# author: Frédéric Darboux <frederic.darboux@orleans.inra.fr> (2012-2015)
# version: 1.06.00
# date: 2015-02-19
#
# License Cecill-V2 <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>
#
# (C) CNRS - Universite d'Orleans - BRGM - INRA (France)
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
SHELL:=/bin/bash

## include specify options for configuration
include ./make_config
###


LIBDIR          := ./lib
BIN             := ./bin
SRC             := $(wildcard *.cpp)
OBJ             := $(SRC:.cpp=.o)

ifeq ($(DEBUG),yes)
	Compil_msg=" Compiling using debug mode "
	Display_command=
else 
	Compil_msg=" Compiling using release mode"
	Display_command= > "/dev/null"
endif

LIBPATH := -L$(LIBDIR)

DIR_SRC :=./Sources

DIR_HEAD :=./Headers

DIR_BENCH := ./Benchmarks

TARGET := FullSWOF_2D

BENCHEXE := BenchEvalFS2D.sh

BENCHPS2DEXE := BenchEvalFS2Dpseudo2D.sh

DIFFCMPEXE := diffCMPFS2D.sh

REPORTDIFFEXE := reportDiffFS2D.sh

LIBRARY := -lschemes -lboundaryconditions -lflux -lfrictions -lrain_infiltration -linitializations -lreconstructions -llimitations -lsave -lparameters -lparser

OBJECTS := FullSWOF_2D.o

INCPATH := -I./$(DIR_HEAD)/libboundaryconditions \
		-I./$(DIR_HEAD)/libflux \
		-I./$(DIR_HEAD)/libfrictions \
		-I./$(DIR_HEAD)/libinitializations \
		-I./$(DIR_HEAD)/liblimitations \
		-I./$(DIR_HEAD)/libparameters \
		-I./$(DIR_HEAD)/libparser \
		-I./$(DIR_HEAD)/librain_infiltration \
		-I./$(DIR_HEAD)/libreconstructions \
		-I./$(DIR_HEAD)/libsave \
		-I./$(DIR_HEAD)/libschemes


###################################################
# Rules to compile all libraries calling Makefiles
###################################################

all: msg construction $(TARGET) clean
	@mv -f $(TARGET) $(BIN)
	@echo " $(TARGET) => $(BIN)"

msg:
	@echo $(Compil_msg)

$(TARGET): $(OBJECTS)
	@$(CPP) -o $@ $^ $(LIBPATH) $(LIBRARY)


%.o: $(DIR_SRC)/%.cpp
	@$(CPP) -o $@ -c $< $(CPPFLAGS) $(INCPATH) 

construction : 
	@cd $(DIR_SRC)/libboundaryconditions && make $(Display_command)
	@cd $(DIR_SRC)/libflux && make $(Display_command)
	@cd $(DIR_SRC)/libfrictions && make $(Display_command)
	@cd $(DIR_SRC)/libinitializations && make $(Display_command)
	@cd $(DIR_SRC)/liblimitations && make $(Display_command)
	@cd $(DIR_SRC)/libparameters && make $(Display_command)
	@cd $(DIR_SRC)/libparser && make $(Display_command)
	@cd $(DIR_SRC)/librain_infiltration && make $(Display_command)
	@cd $(DIR_SRC)/libreconstructions && make $(Display_command)
	@cd $(DIR_SRC)/libsave && make $(Display_command)
	@cd $(DIR_SRC)/libschemes && make $(Display_command)


###################################################
# Rules to remove .o and .a files calling Makefiles
###################################################

clean:
	@cd $(DIR_SRC)/libboundaryconditions && make clean $(Display_command)
	@cd $(DIR_SRC)/libflux && make clean $(Display_command)
	@cd $(DIR_SRC)/libfrictions && make clean $(Display_command)
	@cd $(DIR_SRC)/libinitializations && make clean $(Display_command)
	@cd $(DIR_SRC)/liblimitations && make clean $(Display_command)
	@cd $(DIR_SRC)/libparameters && make clean $(Display_command)
	@cd $(DIR_SRC)/libparser && make clean $(Display_command)
	@cd $(DIR_SRC)/librain_infiltration && make clean $(Display_command)
	@cd $(DIR_SRC)/libreconstructions && make clean $(Display_command)
	@cd $(DIR_SRC)/libsave && make clean $(Display_command)
	@cd $(DIR_SRC)/libschemes && make clean $(Display_command)
	@rm -f $(OBJECTS)

cleanall : clean
	@rm -f ./lib/*.a
	@rm -f $(BIN)/$(TARGET)


###################################################
# Rules to run Benchmarks
###################################################

benchcalc:
	@ # if the code is not compiled, print a message.
	@ if [ ! -x $(BIN)/$(TARGET) ]; then \
		echo "ERROR: $(BIN)/$(TARGET) not found.";\
		echo "Please compile the code by running 'make'. ";	exit 1; \
	fi
	@ echo -e "\n\n***** Emerged bump at rest ******************************"
	@ cd $(DIR_BENCH)/Bump_restemerged && ../../$(BIN)/$(TARGET)
	@ echo -e "\n\n***** MacDonald: Rain fluvial Darcy-Weisbach ************"
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && ../../$(BIN)/$(TARGET)
	@ echo -e "\n\n***** MacDonald: Rain torrential Darcy-Weisbach *********"
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && ../../$(BIN)/$(TARGET)
	@ echo -e "\n\n***** MacDonald: Smooth transition with shock Manning ***"
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && ../../$(BIN)/$(TARGET)
	@ echo -e "\n\n***** MacDonald: Pseudo 2D torrential Manning ***********"
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && ../../$(BIN)/$(TARGET)
	@ echo -e "\n\n***** MacDonald: Pseudo 2D fluvial Manning **************"
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && ../../$(BIN)/$(TARGET)
	@ echo -e "\n\n***** Dry dam break *************************************"
	@ cd $(DIR_BENCH)/Dry_Dam_Break && ../../$(BIN)/$(TARGET)
	@ echo -e "\n\n***** Thacker planar ************************************"
	@ cd $(DIR_BENCH)/Thacker_planar && ../../$(BIN)/$(TARGET)

benchB:
	@ echo " "
	@ echo "WARNING: you will compute all the Benchmark solutions."
	@ echo "This might take a few minutes."
	@ echo " "
	@ sleep 5
	@ # Compute benchmarks
	@ make benchcalc
	@ echo " "
	@ echo "*** Comparison between current and reference results ***"
	@ # Compare results
	@ cd $(DIR_BENCH)/Bump_restemerged && ../../$(BIN)/$(BENCHEXE) analytic.dat 29 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/Bump_restemerged && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/Bump_restemerged && echo -n "Emerged bump at rest:                              " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && ../../$(BIN)/$(BENCHEXE) analytic.dat 33 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && echo -n "MacDonald: Rain fluvial Darcy-Weisbach:            " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && ../../$(BIN)/$(BENCHEXE) analytic.dat 32 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && echo -n "MacDonald: Rain torrential Darcy-Weisbach:         " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && ../../$(BIN)/$(BENCHEXE) analytic.dat 32 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && echo -n "MacDonald: Smooth transition with shock Manning:   " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && ../../$(BIN)/$(BENCHPS2DEXE) analytic.dat 28 Outputs/huz_final.dat 6 9 > comp_USER.dat
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat 
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && echo -n "MacDonald: Pseudo 2D torrential Manning:           " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && ../../$(BIN)/$(BENCHPS2DEXE) analytic.dat 28 Outputs/huz_final.dat 6 99999 > comp_USER.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && echo -n "MacDonald: Pseudo 2D fluvial Manning:              " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/Dry_Dam_Break && ../../$(BIN)/$(BENCHEXE) analytic.dat 26 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/Dry_Dam_Break && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat 
	@ cd $(DIR_BENCH)/Dry_Dam_Break && echo -n "Dry dam break:                                     " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/Thacker_planar && ../../$(BIN)/$(BENCHEXE) analytic.dat 22 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/Thacker_planar && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/Thacker_planar && echo -n "Thacker planar:                                    " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat

benchrefB:
	@ echo " "
	@ echo "REFERENCE solutions for each Benchmark test."
	@ echo "WARNING: the previous reference solutions will be replaced."
	@ echo "This might take a few minutes."
	@ echo "  "
	@ sleep 5
	@ # Compute benchmarks
	@ make benchcalc
	@ echo " "
	@ echo "*** Comparison between reference and standard results ***"
	@ # Get ready for comparisons
	@ cd $(DIR_BENCH) &&\
	for SUBDIR in ./* ; do \
		cd $$SUBDIR && rm -rf Outputs_REFERENCES && mkdir Outputs_REFERENCES && mv Outputs/* Outputs_REFERENCES/ && cd ..; \
	done
	@ # Compare results
	@ cd $(DIR_BENCH)/Bump_restemerged && ../../$(BIN)/$(BENCHEXE) analytic.dat 29 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/Bump_restemerged && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat 
	@ cd $(DIR_BENCH)/Bump_restemerged && echo -n "Emerged bump at rest:                              " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && ../../$(BIN)/$(BENCHEXE) analytic.dat 33 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && echo -n "MacDonald: Rain fluvial Darcy-Weisbach:            " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && ../../$(BIN)/$(BENCHEXE) analytic.dat 32 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && echo -n "MacDonald: Rain torrential Darcy-Weisbach:         " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && ../../$(BIN)/$(BENCHEXE) analytic.dat 32 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && echo -n "MacDonald: Smooth transition with shock Manning:   " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && ../../$(BIN)/$(BENCHPS2DEXE) analytic.dat 28 Outputs_REFERENCES/huz_final.dat 6 9 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat 
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && echo -n "MacDonald: Pseudo 2D torrential Manning:           " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && ../../$(BIN)/$(BENCHPS2DEXE) analytic.dat 28 Outputs_REFERENCES/huz_final.dat 6 99999 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && echo -n "MacDonald: Pseudo 2D fluvial Manning:              " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/Dry_Dam_Break && ../../$(BIN)/$(BENCHEXE) analytic.dat 26 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/Dry_Dam_Break && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/Dry_Dam_Break && echo -n "Dry dam break:                                     " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/Thacker_planar && ../../$(BIN)/$(BENCHEXE) analytic.dat 22 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/Thacker_planar && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/Thacker_planar && echo -n "Thacker planar:                                    " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat

filesfor32bitsbenchmarks:
	@ # get the 32-bit files ready for comparison
	@ echo "*** 32-bit Benchmarks ***" ; \
	cd $(DIR_BENCH) && \
	for SUBDIR in ./* ; do \
		cp $$SUBDIR/comp_STANDARD_32bits.dat $$SUBDIR/comp_STANDARD.dat ; \
	done

filesfor64bitsbenchmarks:
	@ # get the 64-bit files ready for comparison
	@ echo "*** 64-bit benchmarks ***" ; \
	cd $(DIR_BENCH) && \
	for SUBDIR in ./* ; do \
		cp $$SUBDIR/comp_STANDARD_64bits.dat $$SUBDIR/comp_STANDARD.dat ; \
	done

bench32: filesfor32bitsbenchmarks benchB

bench64: filesfor64bitsbenchmarks benchB

benchref32: filesfor32bitsbenchmarks benchrefB

benchref64: filesfor64bitsbenchmarks benchrefB
