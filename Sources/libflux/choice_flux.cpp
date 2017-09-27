/**
 * @file choice_flux.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.01
 * @date 2016-01-04
 *
 * @brief Choice of numerical flux
 * @details 
 * From the value of the corresponding parameter,
 * calls the chosen numerical flux.
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

#include "choice_flux.hpp"


Choice_flux::Choice_flux(int choice){

  /**
   * @details
   * Defines the numerical flux from the value given in the parameters file.
   * @param[in] choice integer that correspond to the chosen numerical flux.
   */

  switch (choice){
    /*  case "Rusanov":*/
    case 1:
      F = new F_Rusanov();
      break;
    /*  case "HLL":*/
    case 2:
      F = new F_HLL();
      break;
    case 3:
      F = new F_HLL2();
      break;
    case 4:
      F = new F_HLLC();
      break;
    case 5:
      F = new F_HLLC2();
      break;
  }
}

void Choice_flux::calcul(SCALAR h_L,SCALAR u_L,SCALAR v_L,SCALAR h_R,SCALAR u_R,SCALAR v_R){
  
  /**
   * @details
   * Calls the calculation of the numerical flux.
   * @param[in] h_L water height at the left of the interface where the flux is calculated.
   * @param[in] u_L velocity (in the x direction) at the left of the interface where the flux is calculated.
   * @param[in] v_L velocity (in the y direction) at the left of the interface where the flux is calculated.
   * @param[in] h_R water height at the right of the interface where the flux is calculated.
   * @param[in] u_R velocity (in the x direction) at the right of the interface where the flux is calculated.
   * @param[in] v_R velocity (in the y direction) at the right of the interface where the flux is calculated.
   */
  
  F->calcul(h_L,u_L,v_L,h_R,u_R,v_R);
}

void Choice_flux::set_tx(SCALAR tx){
  
  /**
   * @details
   * Calls the setting of the value given in parameter to the variable <em> \b tx</em>.
   * @param[in] tx value of dt/dx.
   */
  
  F->set_tx(tx);
}

SCALAR Choice_flux::get_f1(){
  
  /**
   * @details
   * Calls the function to get the first componant of the numerical flux.
   * @return Flux#f1 first componant of the numerical flux.
   */
  
  return F->get_f1();
}

SCALAR Choice_flux::get_f2(){
  
  /**
   * @details
   * Calls the function to get the second componant of the numerical flux.
   * @return Flux#f2 second componant of the numerical flux.
   */
  
  return F->get_f2();
}

SCALAR Choice_flux::get_f3(){
  
  /**
   * @details
   * Calls the function to get the third componant of the numerical flux.
   * @return Flux#f3 third componant of the numerical flux.
   */
  
  return F->get_f3();
}

SCALAR Choice_flux::get_cfl(){
  
  /**
   * @details
   * Calls the function to get the value of the CFL.
   * @return Flux#cfl value of the CFL.
   */
  
  return F->get_cfl();
}

Choice_flux::~Choice_flux(){
  if (F!=NULL){
    delete F;
    F = NULL;
  }
}
