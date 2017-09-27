/**
 * @file bc_neumann.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Neumann condition
 * @details 
 * Boundary condition:
 * Neumann condition (the normal derivative is null).
 *
 * @copyright License Cecill-V2 \n
 * <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>
 *
 * (c) CNRS - Universite d'Orleans - BRGM (France)
 */
/* 
 *
 * This file is part of FullSWOF_2D software. 
 * <https://sourcesup.renater.fr/projects/fullswof-2d/> 
 *
 * FullSWOF_2D = Full Shallow-Water equations for Overland Flow, 
 * in two dimensions of space.
 * This software is a computer program whose purpose is to compute
 * solutions for 2D Shallow-Water equations.
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


#include "bc_neumann.hpp"

Bc_Neumann::Bc_Neumann(Parameters & par,TAB & z,int n1, int n2): Boundary_condition(par){
  
  /**
   * @details
   * @param[in] par parameter, contains all the values from the parameters file (unused).
   * @param[in] n1 integer to specify whether it is the left (-1) or the right (1) boundary.
   * @param[in] n2 integer to specify whether it is the bottom (-1) or the top (1) boundary.
   * @param[in,out] z vector that represents the topography with suitable values on the fictive cells. 
   */
  
  unormbound = 0.;
  utanbound = 0.;
  hbound = 0.;
  
  // Initialization of the topography at the boundary
  if ((-1==n1)&&(0==n2)){
    // Left Boundary condition
    for (int j=1 ; j<NYCELL+1 ; j++){
      //      z[0][j]=z[1][j];
      z[0][j]=z[1][j]-((z[NXCELL][j]-z[1][j])/(NXCELL-1));
    }
  }
  if ((0==n1)&&(-1==n2)){
    // Bottom Boundary condition
    for (int i=1 ; i<NXCELL+1 ; i++){
      //      z[i][0]=z[i][1];
      z[i][0]=z[i][1]-((z[i][NYCELL]-z[i][1])/(NYCELL-1));
    }
  }
  if ((1==n1)&&(0==n2)) {
    // Right Boundary condition
    for (int j=1 ; j<NYCELL+1 ; j++){
      //      z[NXCELL+1][j]=z[NXCELL][j];
      z[NXCELL+1][j]=z[NXCELL][j]+((z[NXCELL][j]-z[1][j])/(NXCELL-1));
    }
  }
  if ((0==n1)&&(1==n2)) {
    // Top Boundary condition
    for (int i=1 ; i<NXCELL+1 ; i++){
      //      z[i][NYCELL+1]=z[i][NYCELL];
      z[i][NYCELL+1]=z[i][NYCELL]+((z[i][NYCELL]-z[i][1])/(NYCELL-1));
    }
  }
}

void Bc_Neumann::calcul(SCALAR hin,SCALAR unorm_in,SCALAR utan_in,SCALAR hfix,SCALAR qfix,SCALAR hin_oppbound,SCALAR unorm_in_oppbound,SCALAR utan_in_oppbound,SCALAR time, int n1,int n2){
  
  /**
   * @details
   * @param[in] hin water height of the first cell inside the domain.
   * @param[in] unorm_in normal velocity of the first cell inside the domain.
   * @param[in] utan_in tangential velocity of the first cell inside the domain.
   * @param[in] hfix fixed (imposed) value of the water height (unused).
   * @param[in] qfix fixed (imposed) value of the discharge (unused).
   * @param[in] hin_oppbound value of the water height of the first cell inside the domain at the opposite bound (unused).
   * @param[in] unorm_in_oppbound value of the normal velocity of the first cell inside the domain at the opposite bound (unused).
   * @param[in] utan_in_oppbound value of the tangential velocity of the first cell inside the domain at the opposite bound (unused).
   * @param[in] time current time (unused).
   * @param[in] n1 integer to specify whether it is the left (-1) or the right (1) boundary (unused).
   * @param[in] n2 integer to specify whether it is the bottom (-1) or the top (1) boundary (unused).
   * @par Modifies
   * Boundary_condition#hbound water height on the fictive cell.\n
   * Boundary_condition#unormbound normal velocity on the fictive cell.\n
     * Boundary_condition#utanbound tangential velocity on the fictive cell.\n
   */
  
  (void) qfix; //unused variable
  (void) hfix; //unused variable
  (void) hin_oppbound; //unused variable
  (void) unorm_in_oppbound; //unused variable
  (void) utan_in_oppbound; //unused variable
  (void) time; //unused variable
  (void) n1; //unused variable
  (void) n2; //unused variable
  
  unormbound = unorm_in;
  utanbound = utan_in;
  hbound = hin;
}

Bc_Neumann::~Bc_Neumann(){
}
