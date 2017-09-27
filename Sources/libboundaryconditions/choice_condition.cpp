/**
 * @file choice_condition.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Choice of boundary condition
 * @details 
 * From the value of the corresponding parameter,
 * calls the chosen boundary condition.
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


#include "choice_condition.hpp"

Choice_condition::Choice_condition(int choice,Parameters & par, TAB & z,int n1, int n2){
  
  /**
   * @details
   * Defines the boundary condition from the value given in the parameters file.
   * @param[in] choice integer that correspond to the chosen boundary condition.
   * @param[in] par parameter, contains all the values from the parameters file.
   * @param[in] z array that represents the topography.
   * @param[in] n1 integer to specify whether it is the left (-1) or the right (1) boundary.
   * @param[in] n2 integer to specify whether it is the bottom (-1) or the top (1) boundary.
   */
  
  switch (choice){
    case 1:
      C = new Bc_imp_height(par,z,n1,n2);
      break;
    case 2:
      C = new Bc_wall(par,z,n1,n2);
      break;
    case 3:
      C = new Bc_Neumann(par,z,n1,n2);
      break;
    case 4:
      C = new Bc_periodic(par,z,n1,n2);
      break;
    case 5:
      C = new Bc_imp_discharge(par,z,n1,n2);
      break;
  }
}

void Choice_condition::calcul(SCALAR hin, SCALAR unorm_in, SCALAR utan_in, SCALAR hfix, SCALAR qfix, SCALAR hin_oppbound, SCALAR unorm_in_oppbound, SCALAR utan_in_oppbound, SCALAR time, int n1, int n2){
  
  /**
   * @details
   * Calls the calculation of the boundary condition.
   * @param[in] hin water height of the first cell inside the domain.
   * @param[in] unorm_in normal velocity of the first cell inside the domain.
   * @param[in] utan_in tangential velocity of the first cell inside the domain.
   * @param[in] hfix fixed (imposed) value of the water height.
   * @param[in] qfix fixed (imposed) value of the discharge.
   * @param[in] hin_oppbound value of the water height of the first cell inside the domain at the opposite bound.
   * @param[in] unorm_in_oppbound value of the normal velocity of the first cell inside the domain at the opposite bound.
   * @param[in] utan_in_oppbound value of the tangential velocity of the first cell inside the domain at the opposite bound.
   * @param[in] time current time.
   * @param[in] n1 integer to specify whether it is the left (-1) or the right (1) boundary.
   * @param[in] n2 integer to specify whether it is the bottom (-1) or the top (1) boundary.
   */
  
  C->calcul(hin,unorm_in,utan_in,hfix,qfix,hin_oppbound,unorm_in_oppbound,utan_in_oppbound,time, n1, n2);
}

SCALAR Choice_condition::get_hbound(){
  
  /**
   * @details
   * Calls the function to get the water height on the fictive cell.
   * @return Boundary_condition#hbound water height on the fictive cell for the chosen boundary condition.
   */
  
  return C->get_hbound();
}

SCALAR Choice_condition::get_unormbound(){
  
  /**
   * @details
   * Calls the function to get the normal velocity on the fictive cell.
   * @return Boundary_condition#unormbound normal velocity on the fictive cell for the chosen boundary condition.
   */
  
  return C->get_unormbound();
}

SCALAR Choice_condition::get_utanbound(){
  
  /**
   * @details
   * Calls the function to get the tangential velocity on the fictive cell.
   * @return Boundary_condition#utanbound tangential velocity on the fictive cell for the chosen boundary condition.
   */
  
  return C->get_utanbound();
}

Choice_condition::~Choice_condition()
{
  if (C!=NULL) 
  {
    delete C;
    C=NULL;
  }
}
