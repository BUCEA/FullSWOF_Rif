/**
 * @file choice_init_huv.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Choice of initialization for h, u and v
 * @details 
 * From the value of the corresponding parameter,
 * calls the chosen initialization of the water height and of the velocity.
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

#include "choice_init_huv.hpp"

Choice_init_huv::Choice_init_huv(Parameters & par){
  
  /**
   * @details
   * Defines the initialization of the water height and of the velocity
   * from the value given in the parameters file.
   * @param[in] par parameter, contains all the values from the parameters file.
   */
  
  switch (par.get_huv()){
    case 1:
      huv_init = new Huv_read(par);
      break;
    case 2:
      huv_init = new Huv_generated(par);
      break;
    case 3:
      huv_init = new Huv_generated_Thacker(par);
      break;
    case 4:
      huv_init = new Huv_generated_Radial_Dam_dry(par);
      break;
    case 5:
      huv_init = new Huv_generated_Radial_Dam_wet(par);
      break;
  }
}

void Choice_init_huv::initialization(TAB & h,TAB & u,TAB & v){
  
  /**
   * @details
   * Calls the initialization of the water height and of the velocity.
   * @param[in] h water height.
   * @param[in] u first componant of the velocity.
   * @param[in] v second componant of the velocity.
   */
  
  huv_init->initialization(h,u,v);
}

Choice_init_huv::~Choice_init_huv(){
  delete huv_init;
}

