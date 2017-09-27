/**
 * @file no_infiltration.cpp
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief No infiltration
 * @details
 * %Infiltration: there is no infiltration.
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

#include "no_infiltration.hpp"

No_Infiltration::No_Infiltration(Parameters & par):Infiltration(par){
  
  /**
   * @details
   * @param[in] par parameter, contains all the values from the parameters file (unused).
   */
}


void No_Infiltration:: calcul(const TAB & h,const TAB & Vin_tot,SCALAR dt){
  
  /**
   * @details No computation (water height and infiltrated volume remain unchanged).
     * @param[in] h water height.
     * @param[in] Vin_tot total infiltrated volume.
     * @param[in] dt time step (unused).
     * @par Modifies
     * Infiltration#hmod modified water height.\n
     * Infiltration#Vin total infiltrated volume containing the current time step.
     */
  
  (void) dt; //unused variable
  
  for (int i=1 ; i<NXCELL+1 ; i++){
    for (int j=1 ; j<NYCELL+1 ; j++){
      hmod[i][j] = h[i][j]; // water height is no modified.
      Vin[i][j] = Vin_tot[i][j]; // total infiltrated volum is no modified.     
    } //end for j
  } //end for i
}


No_Infiltration::~No_Infiltration(){
}



