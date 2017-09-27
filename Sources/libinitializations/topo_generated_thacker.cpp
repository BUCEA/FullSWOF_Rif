/**
 * @file topo_generated_thacker.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Thacker configuration
 * @details 
 * Initialization of the topography:
 * topography with a shape of a paraboloid of revolution for Thacker's Benchmark.
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

#include "topo_generated_thacker.hpp"

Topo_generated_Thacker::Topo_generated_Thacker(Parameters & par):Initialization_topo(par){
  
  /**
   * @details   
   * Defines the parameters of the paraboloid.
   * @param[in] par parameter, contains all the values from the parameters file (unused).
   */
  
  a=1.;
  z0=0.1;
  r0=0.8;
  Aa=((a*a)-(r0*r0))/((a*a)+(r0*r0));
  omega=sqrt(8*GRAV*z0)/a;
  mid_x = (NXCELL*DX)/2.;
  mid_y = (NYCELL*DX)/2.;
}


void Topo_generated_Thacker::initialization(TAB & topo){
  
  /**
   * @details
   * Initializes the topography to \f$ h_0\left(\frac{(x-Lx/2)^2+(y-Ly/2)^2}{a^2}-1\right) \f$, see \cite Thacker81.
   * @param[in, out] topo topography.
   */
  
  for (int i=1 ; i<NXCELL+1 ; i++){
    for (int j=1 ; j<NYCELL+1 ; j++){
      x=(i-0.5)*DX;
      y=(j-0.5)*DY;
      
      radius2=(x-mid_x)*(x-mid_x)+(y-mid_y)*(y-mid_y);
      
      topo[i][j] = -z0*(1.-(radius2)/(a*a));
    }
  }
}

Topo_generated_Thacker::~Topo_generated_Thacker(){
}
