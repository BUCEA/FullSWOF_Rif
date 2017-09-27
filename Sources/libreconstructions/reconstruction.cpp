/**
 * @file reconstruction.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief %Reconstruction
 * @details 
 * Common part for all the reconstructions.
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

#include "reconstruction.hpp"

Reconstruction::Reconstruction(Parameters & par, TAB & z):NXCELL(par.get_Nxcell()),NYCELL(par.get_Nycell()){
  
  /**
   * @details
   * Defines the number of cells, the slope limiter, and initializes Reconstruction#z1l, Reconstruction#z2l, Reconstruction#z1r, Reconstruction#z2r, Reconstruction#delta_z1, Reconstruction#delta_z2.
   * @param[in] par parameter, contains all the values from the parameters file.
   * @param[in] z topography.
   */
  
  limiter = new Choice_limiter(par.get_lim());
  
  z1r.resize(NXCELL+1); // i : 0->NXCELL
  z1l.resize(NXCELL+2); // i : 1->NXCELL+1
  z2r.resize(NXCELL+1); // i : 1->NXCELL
  z2l.resize(NXCELL+1); // i : 1->NXCELL
  delta_z1.resize(NXCELL+1); // i : 0->NXCELL
  delta_z2.resize(NXCELL+1); // i : 1->NXCELL
  
  z1r[0].resize(NYCELL+1); // lj: 1->NYCELL
  delta_z1[0].resize(NYCELL+1); // j : 1->NYCELL
  
  
  for (int i=1 ; i<NXCELL+1 ; i++){
    z1r[i].resize(NYCELL+1); // j : 1->NYCELL
    z1l[i].resize(NYCELL+1); // j : 1->NYCELL
    z2r[i].resize(NYCELL+1); // j : 0->NYCELL
    z2l[i].resize(NYCELL+2); // j : 1->NYCELL+1
    delta_z1[i].resize(NYCELL+1); // j : 1->NYCELL
    delta_z2[i].resize(NYCELL+1); // j : 0->NYCELL
  }
  
  z1l[NXCELL+1].resize(NYCELL+1);
  
  for (int j=1 ; j<NYCELL+1 ; j++){
    z1r[0][j] = z[0][j];
    z1l[NXCELL+1][j] = z[NXCELL+1][j];
    delta_z1[0][j] = z[1][j]-z[0][j];
  } //end for j
  
  for (int i=1 ; i<NXCELL+1 ; i++){
    z2r[i][0] = z[i][0];
    z2l[i][NYCELL+1] = z[i][NYCELL+1];
    delta_z2[i][0] = z[i][1]-z[i][0];
    for (int j=1 ; j<NYCELL+1 ; j++){
      delta_z1[i][j] = z[i+1][j]-z[i][j];
      delta_z2[i][j] = z[i][j+1]-z[i][j];
    } //end for j
  } //end for i
}

Reconstruction::~Reconstruction(){
  
  if (limiter != NULL){
    delete limiter;
    limiter = NULL;
  }
  
  z1r[0].clear();
  delta_z1[0].clear();
  
  for (int i=1 ; i<=NXCELL ; i++){
    z1r[i].clear();
    z1l[i].clear();
    z2r[i].clear();
    z2l[i].clear();
    delta_z1[i].clear();
    delta_z2[i].clear();
  }
  
  z1l[NXCELL+1].clear();
  
  z1r.clear();
  z1l.clear();
  z2r.clear();
  z2l.clear();
  delta_z1.clear();
  delta_z2.clear();
}

