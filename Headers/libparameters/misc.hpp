/**
 * @file misc.hpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.07.01            //Revision 358, 06/01/2016
 * @date 2016-06-01            //Revision 358, 06/01/2016
 *
 * @brief Definitions
 * @details 
 * Defines the constants, the types used in the code and contains the `include'.
 *
 * @copyright License Cecill-V2  \n
 * <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>
 *
 * (c) CNRS - Universite d'Orleans - INRA (France)
 */
/*
 * This file is part of FullSWOF_2D software. 
 * <https://sourcesup.renater.fr/projects/fullswof-1d/> 
 *
 * FullSWOF_2D = Full Shallow-Water equations for Overland Flow, 
 * in one dimension of space.
 * This software is a computer program whose purpose is to compute
 * solutions for 1D Shallow-Water equations.
 *
 * LICENSE
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * <http://www.cecill.info>.
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading, using, modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean that it is complicated to manipulate, and that also
 * therefore means that it is reserved for developers and experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and, more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *
 ******************************************************************************/

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <cfloat>
#include <unistd.h>
#include <ctime>

#define max(a,b) (a>=b?a:b)
#define min(a,b) (a<=b?a:b)

#define GRAV 9.81
#define GRAV_DEM 4.905
#define CONST_CFL_X 0.5
#define CONST_CFL_Y  0.5
#define HE_CA 1.e-12
#define VE_CA 1.e-12
#define MAX_CFL_X 0.
#define MAX_CFL_Y 0.
#define MAX_ITER 1000000000
/*NB_CHAR is the maximum length of a comment line */
#define NB_CHAR 256
#define ZERO 0.
#define IE_CA 1.e-8
#define EPSILON 1.e-13
#define VERSION "FullSWOF_2D version 1.07.01, 2016-06-01"         //Revision 358, 06/01/2016

//RATIO_CLOSE_CELL is used to verify that the input data is very close to cell center
#define RATIO_CLOSE_CELL 1.e-3

using namespace std;

typedef double SCALAR;/*If you change the data type of SCALAR, don't forget to change the macro definition MAX_SCAL*/

/*Maximum finite representable floating-point number of SCALAR.*/
#define MAX_SCAL DBL_MAX //MAX_SCAL=-1
//predefine two dimensional vector constructed from scalar s (fills all components with s) 
typedef vector< vector< SCALAR > > TAB;
//predefine one dimensional vector constructed from scalar s (fills all components with s) 
typedef vector<SCALAR> vect;