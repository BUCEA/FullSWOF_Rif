/** 
 * @file no_friction.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief No friction
 * @details 
 * @brief Friction law: does no computation.
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

#include "no_friction.hpp"

No_Friction::No_Friction(Parameters & par):Friction(par){
  
  /**
   * @details
   * @param[in] par parameter, contains all the values from the parameters file.
   */
  
}

void No_Friction::calcul(const TAB & uold,const TAB &  vold,const TAB &  hnew, const TAB & q1new,const TAB &  q2new, SCALAR dt){
  
  /**
   * @details
   * No computation (no friction).
   * @param[in] uold velocity in the first direction at the previous time (\f$n\f$ if you are calculating the \f$n+1\f$th time step) (unused).
   * @param[in] vold velocity in the second direction at the previous time (\f$n\f$ if you are calculating the \f$n+1\f$th time step) (unused).
   * @param[in] hnew water height after the Shallow-Water computation (without friction) (unused).
   * @param[in] q1new discharge in the first direction after the Shallow-Water computation (without friction) (unused).
   * @param[in] q2new discharge in the second direction after the Shallow-Water computation (without friction) (unused).
   * @param[in] dt time step (unused).
   * @par Modifies
   * Friction#q1mod discharge in the first direction modified by the friction term,\n
   * Friction#q2mod discharge in the second direction modified by the friction term.
   */
  
  (void) uold; //unused variable
  (void) vold; //unused variable
  (void) hnew; //unused variable
  (void) dt; //unused variable
  
  for (int i=1 ; i<NXCELL+1 ; i++){
    for (int j=1 ; j<NYCELL+1 ; j++){
      q1mod[i][j] = q1new[i][j];
      q2mod[i][j] = q2new[i][j];
    } //end for j
  } //end for i
}

void No_Friction::calculSf(const TAB & h,const TAB &  u,const TAB &  v){
  
  /**
   * @details
   * Explicit friction term: \f$ S_f=0\f$.
   * @param[in] h water height (unused).
   * @param[in] u velocity in the first direction (unused).
   * @param[in] v velocity in the second direction (unused).
   * @par Modifies
   * Friction#Sf1 explicit friction term in the first direction,\n
   * Friction#Sf2 explicit friction term in the second direction.
   * @note This explicit friction term will be used for erosion.
   */
  
  (void) h; //unused variable
  (void) u; //unused variable
  (void) v; //unused variable
  
  for (int i=1 ; i<NXCELL+1 ; i++){
    for (int j=1 ; j<NYCELL+1 ; j++){
      Sf1[i][j] = ZERO;
      Sf2[i][j] = ZERO;
    } //end for j
  } //end for i
}

No_Friction::~No_Friction(){
}

