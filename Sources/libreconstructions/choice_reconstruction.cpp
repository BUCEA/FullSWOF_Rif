/**
 * @file choice_reconstruction.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Choice of reconstruction
 * @details 
 * From the value of the corresponding parameter,
 * calls the chosen reconstruction.
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

#include "choice_reconstruction.hpp"

Choice_reconstruction::Choice_reconstruction(Parameters & par, TAB & z){
  
  /**
   * @details
   * Defines the reconstruction from the value given in the parameters file.
   * @param[in] par parameter, contains all the values from the parameters file.
   * @param[in] z array that represents the topography.
   */
  
  switch (par.get_rec()){
    case 1:
      Rec = new MUSCL(par,z);
      break;
    case 2:
      Rec = new ENO(par,z);
      break;
    case 3:
      Rec = new ENO_mod(par,z);
      break;
    }
}

void Choice_reconstruction::calcul(TAB & h,TAB & u,TAB & v,TAB & z,TAB & delzc1,TAB & delzc2,TAB & delz1,TAB & delz2,TAB & h1r,TAB & u1r,TAB & v1r,TAB & h1l,TAB & u1l,TAB & v1l,TAB & h2r,TAB & u2r,TAB & v2r,TAB & h2l,TAB & u2l,TAB & v2l){
  
  /**
   * @details
   * Calls the calculation of the second order reconstruction in space.
   * @param[in] h water height.
   * @param[in] u velocity of the flow in the first direction.
   * @param[in] v velocity of the flow in the second direction.
   * @param[in] z topography.
   * @param[out] delzc1 difference between the reconstructed topographies on the left and on the right boundary of a cell in the first direction.
   * @param[out] delzc2 difference between the reconstructed topographies on the left and on the right boundary of a cell in the second direction.
   * @param[out] delz1 difference between two reconstructed topographies on the same boundary (from two adjacent cells) in the first direction.
   * @param[out] delz2 difference between two reconstructed topographies on the same boundary (from two adjacent cells) in the seond direction.
   * @param[out] h1r reconstructed water height on the right of the cell in the first direction.
   * @param[out] u1r first componant of the reconstructed velocity on the right of the cell in the first direction.
   * @param[out] v1r second componant of the reconstructed velocity on the right of the cell in the first direction.
   * @param[out] h1l reconstructed water height on the left of the cell in the first direction.
   * @param[out] u1l first componant of the reconstructed velocity on the left of the cell in the first direction.
   * @param[out] v1l second componant of the reconstructed velocity on the left of the cell in the first direction.
   * @param[out] h2r reconstructed water height on the right of the cell in the second direction.
   * @param[out] u2r first componant of the reconstructed velocity on the right of the cell in the second direction.
   * @param[out] v2r second componant of the reconstructed velocity on the right of the cell in the second direction.
   * @param[out] h2l reconstructed water height on the left of the cell in the second direction.
   * @param[out] u2l first componant of the reconstructed velocity on the left of the cell in the second direction.
   * @param[out] v2l second componant of the reconstructed velocity on the left of the cell in the second direction.
   */
  
  Rec->calcul(h,u,v,z,delzc1,delzc2,delz1,delz2,h1r,u1r,v1r,h1l,u1l,v1l,h2r,u2r,v2r,h2l,u2l,v2l);
}

Choice_reconstruction::~Choice_reconstruction(){
  if (Rec != NULL){
    delete Rec;
    Rec = NULL;
  }
}
