/**
 * @file choice_init_rif.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Choice of initialization for r, i and f
 * @details 
 * From the value of the corresponding parameter,
 * calls the chosen initialization of the rainfall, inifiltration and friction choice.
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

#include "choice_init_rif.hpp"

Choice_init_rif::Choice_init_rif(Parameters & par)
{
  
  /**
   * @details
   * Defines the initialization of the rainfall, infiltration and friction choice
   * from the value given in the parameters file.
   * @param[in] par parameter, contains all the values from the parameters file.
   */
  
  switch (par.get_rif())
  {
    case 1:
      rif_init = new Rif_read(par);
      break;
    /*case 2:
      rif_init = new Rif_generated(par);
      break;
    case 3:
      rif_init = new Rif_generated_Thacker(par);
      break;
    case 4:
      rif_init = new Rif_generated_Radial_Dam_dry(par);
      break;
    case 5:
      rif_init = new Rif_generated_Radial_Dam_wet(par);
      break;*/
  }
}

void Choice_init_rif::initialization(TAB & rain_c,TAB & infi_c,TAB & fric_c){
  
  /**
   * @details
   * Calls the initialization of the rainfall, infiltration and friction choice.
   * @param[in] rain_c rainfall choice.
   * @param[in] infi_c infiltration choice.
   * @param[in] fric_c friction choice.
   */
  
  rif_init->initialization(rain_c,infi_c,fric_c);
}

Choice_init_rif::~Choice_init_rif(){
  delete rif_init;
}

