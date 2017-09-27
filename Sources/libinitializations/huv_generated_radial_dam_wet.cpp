/**
 * @file huv_generated_radial_dam_wet.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Wet radial dam break configuration
 * @details 
 * Initialization of the water height and of the velocity:
 * case of a radial dam break on a wet domain.
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

#include "huv_generated_radial_dam_wet.hpp"

Huv_generated_Radial_Dam_wet::Huv_generated_Radial_Dam_wet(Parameters & par):Initialization_huv(par){
  
  /**
   * @details   
   * Defines the position of the dam (half of the domain), the water height before the dam (5 millimeters),
   * the water height after the dam (4 millimeter) and the velocity (0 m/s), see \cite Goutal97, \cite Audusse00. 
   * @param[in] par parameter, contains all the values from the parameters file (unused).
   */
  
  dvalu = dvalv = 0.;
  dvalh = 0.001;
  dvalhh=0.004;
  
  nmid_x = NXCELL/ 2.;
  nmid_y = NYCELL/ 2.;
  
  nradius = 0.1*NXCELL;
}

void Huv_generated_Radial_Dam_wet::initialization(TAB & h,TAB & u,TAB & v){
  
  /**
   * @details
   * Initializes the water height and the velocity, before and after the dam.
   * @param[in, out] h water height.
   * @param[in, out] u first componant of the velocity.
   * @param[in, out] v second componant of the velocity.
   */
  
  for (int i=0 ; i<=NXCELL+1 ; i++){
    for (int j=0 ; j<=NYCELL+1 ; j++){
      h[i][j]=dvalh;
      u[i][j]=dvalu;
      v[i][j]=dvalv;
      if ((i-nmid_x)*(i-nmid_x)+(j-nmid_y)*(j-nmid_y) <= nradius*nradius)
        h[i][j] = dvalh+dvalhh;
    }
  }
}

Huv_generated_Radial_Dam_wet::~Huv_generated_Radial_Dam_wet(){
}
