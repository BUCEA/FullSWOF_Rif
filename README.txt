#### FullSWOF_2D
#### version 1.07.00, 2016-03-14

FullSWOF_2D stands for ``Full Shallow Water equations for Overland Flow in two 
dimension of space''. In this software, the Shallow Water (or Saint-Venant) equations are solved 
using finite volumes and numerical methods especially chosen for hydrodynamic purposes 
(transitions between wet and dry areas, small water heights, steady preservation ...). 

This software is distributed under CeCILL-V2 (GPL compatible) free software license
(see <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>).

The following explains how to compile FullSWOF_2D.
A more complete documentation can be found in the doc/Documentation.pdf file.

##############################
Users
##############################

Windows' users: look at the application note 
<https://sourcesup.renater.fr/docman/view.php/895/3949/AppNote-windows.pdf>

When you are in the FullSWOF_2D directory, write the following lines:

make cleanall
make
cd Exp01
../bin/FullSWOF_2D

All the results are saved in the Exp01/Outputs/ directory, in .dat files.
Your "parameters.txt" file and your possible topography/huv files must be in the Exp01/Inputs/ directory.

Then check for proper functionning by running once
make benchref32
	or
make benchref64
(depending if your operating system is 32 or 64 bits).

See the documentation for more details.

##############################
Developers
##############################

You should start using FullSWOF_2D by running all the benchmarks and saving the reference solutions.
This will allow you to compare the present results with your future results after you modify the source files.
To do so, run once
make benchref32
	or
make benchref64
(depending if your operating system is 32 or 64 bits).
After, each time you want to check if your changes in the source files affect the results, run
make bench32
	or
make bench64
(depending if your operating system is 32 or 64 bits).

See the documentation for more details.

Always comment the files, at the beginning of the file, using Doxygen syntax (www.doxygen.org/).

---------- Doxygen -----------

Note: To simplify these operations, you can run the script UpdateDateVersion located in the bin folder.

In order to generate the Doxygen html file, the Doxygen_config_file_html file is saved in the doc directory.
To run Doxygen, from the FullSWOF_2D directory, use the command:

doxygen doc/Doxygen_config_file_html

Warning: Graphviz (http://www.graphviz.org/) must be in your PATH to generate HTML diagrams. 
If not, change the HAVE_DOT parameter of the Doxygen_config_file_html.

=> In the doc/html/ directory, index.html is created.

To generate the Doxygen latex (pdf) file, you must use the Doxygen_config_file_latex file and compile the tex file:

doxygen doc/Doxygen_config_file_latex
cd doc/latex
make

=> In the doc/latex directory, refman.pdf is created.

-------- Check list  --------

Before creating a new tag / version, check the compilation for errors and warning messages. Then run the benchmarks under various operating systems.
Also detail the modifications in the changelog file.
Then, please refer to the Doc-NewReleaseHowTo file in the doc folder to create a new version and a new tag.
	
