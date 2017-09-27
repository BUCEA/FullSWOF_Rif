/**
 * @file eno_mod.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Modified %ENO reconstruction
 * @details 
 * Linear reconstruction: modified %ENO.
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

#include "eno_mod.hpp"


ENO_mod::ENO_mod(Parameters & par,TAB & z):Reconstruction(par, z),AMORTENO(par.get_amortENO()),MODIFENO(par.get_modifENO()){
  
  /**
   * @details
   * Initializations.
   * @param[in] par parameter, contains all the values from the parameters file.
   * @param[in] z topography.
   */
  
  som_z1.resize(NXCELL+1); // i : 1->NXCELL
  som_z2.resize(NXCELL+1); // i : 1->NXCELL
  
  for (int i=1 ; i<NXCELL+1 ; i++){
    som_z1[i].resize(NYCELL+1); // j : 1->NYCELL
    som_z2[i].resize(NYCELL+1); // j : 1->NYCELL
  }
  
  for (int i=1 ; i<NXCELL+1 ; i++){
    for (int j=1 ; j<=NYCELL ; j++){
      som_z1[i][j] = z[i-1][j]-2.*z[i][j]+z[i+1][j];
      som_z2[i][j] = z[i][j-1]-2.*z[i][j]+z[i][j+1];
    } //end for j
  } //end for i
}

void ENO_mod::calcul(TAB & h,TAB & u,TAB & v,TAB & z,TAB & delzc1,TAB & delzc2,TAB & delz1,TAB & delz2,TAB & h1r,TAB & u1r,TAB & v1r,TAB & h1l,TAB & u1l,TAB & v1l,TAB & h2r,TAB & u2r,TAB & v2r,TAB & h2l,TAB & u2l,TAB & v2l){
  
  /**
   * @details
   * Calls the calculation of the second order reconstruction in space, with a modified %ENO formulation, see \cite Bouchut04, \cite Bouchut07b.
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
  
  
  for (int j=1 ; j<NYCELL+1 ; j++){
    
    ddh1=0.;
    ddz1=0.;
    ddu1=0.;
    ddv1=0.;
    
    hh1 = h[0][j]-2.*h[1][j]+h[2][j];
    uu1 = u[0][j]-2.*u[1][j]+u[2][j];
    vv1 = v[0][j]-2.*v[1][j]+v[2][j];
    
    delta_h1 = h[1][j]-h[0][j];
    delta_u1 = u[1][j]-u[0][j];
    delta_v1 = v[1][j]-v[0][j];
    
    for (int i=1 ; i<NXCELL ; i++){
      
      hh2 = h[i][j]-2.*h[i+1][j]+h[i+2][j];
      uu2 = u[i][j]-2.*u[i+1][j]+u[i+2][j];
      vv2 = v[i][j]-2.*v[i+1][j]+v[i+2][j];
      
      limiter->calcul(hh1,hh2);
      ddh2 = AMORTENO*limiter->get_rec();
      limiter->calcul(uu1,uu2);
      ddu2 = AMORTENO*limiter->get_rec();
      limiter->calcul(hh1+som_z1[i][j],hh2+som_z1[i+1][j]);
      ddz2 = AMORTENO*limiter->get_rec();
      limiter->calcul(vv1,vv2);
      ddv2 = AMORTENO*limiter->get_rec();
      
      hh1 = hh2;
      uu1 = uu2;
      vv1 = vv2;
      
      delta_h2 = h[i+1][j]-h[i][j];
      delta_u2 = u[i+1][j]-u[i][j];
      delta_v2 = v[i+1][j]-v[i][j];
      
      limiter->calcul(delta_h1+ddh1*0.5,delta_h2-ddh2*0.5);
      a2 = limiter->get_rec(); 
      limiter->calcul(delta_h1+delta_z1[i-1][j]+ddz1*0.5,delta_h2+delta_z1[i][j]-ddz2*0.5);
      a4 = limiter->get_rec(); 
      
      limiter->calcul(delta_h1,delta_h2);
      a1 = limiter->get_rec();
      limiter->calcul(2*MODIFENO*a1,a2);
      dh = limiter->get_rec();
      
      limiter->calcul(delta_h1+delta_z1[i-1][j],delta_h2+delta_z1[i][j]);
      a3 = limiter->get_rec();
      limiter->calcul(2*MODIFENO*a3,a4);
      dz_h = limiter->get_rec();
      
      limiter->calcul(delta_u1+ddu1*0.5,delta_u2-ddu2*0.5);
      du = limiter->get_rec();
      limiter->calcul(delta_v1+ddv1*0.5,delta_v2-ddv2*0.5);
      dv = limiter->get_rec();
      
      delta_h1 = delta_h2;
      delta_u1 = delta_u2;
      delta_v1 = delta_v2;
      
      h1r[i][j]=h[i][j]+dh*0.5;
      h1l[i][j]=h[i][j]-dh*0.5;
      
      z1r[i][j]=z[i][j]+0.5*(dz_h-dh);
      z1l[i][j]=z[i][j]+0.5*(dh-dz_h);
      delzc1[i][j] = z1r[i][j]-z1l[i][j];
      delz1[i-1][j] = z1l[i][j]-z1r[i-1][j];
      
      if (h[i][j]>0.){
        u1r[i][j]=u[i][j]+h1l[i][j]*du*0.5/h[i][j];
        u1l[i][j]=u[i][j]-h1r[i][j]*du*0.5/h[i][j];
        v1r[i][j]=v[i][j]+h1l[i][j]*dv*0.5/h[i][j];
        v1l[i][j]=v[i][j]-h1r[i][j]*dv*0.5/h[i][j];
      }
      else{
        u1r[i][j]=u[i][j]+du*0.5;
        u1l[i][j]=u[i][j]-du*0.5;
        v1r[i][j]=v[i][j]+dv*0.5;
        v1l[i][j]=v[i][j]-dv*0.5;
      } //end if
      
      ddh1=ddh2;
      ddz1=ddz2;
      ddu1=ddu2;
      ddv1=ddv2;
      
    } //end for
    
    
    limiter->calcul(delta_h1+ddh1*0.5,h[NXCELL+1][j]-h[NXCELL][j]);
    a2 = limiter->get_rec();
    limiter->calcul(delta_h1+delta_z1[NXCELL-1][j]+ddz1*0.5,h[NXCELL+1][j]-h[NXCELL][j]+delta_z1[NXCELL][j]);
    a4 = limiter->get_rec();
    
    limiter->calcul(delta_h1,h[NXCELL+1][j]-h[NXCELL][j]);
    a1 = limiter->get_rec();
    limiter->calcul(2*MODIFENO*a1,a2);
    dh = limiter->get_rec();
    
    limiter->calcul(delta_h1+delta_z1[NXCELL-1][j],h[NXCELL+1][j]-h[NXCELL][j]+delta_z1[NXCELL][j]);
    a3 = limiter->get_rec();
    limiter->calcul(2*MODIFENO*a3,a4);
    dz_h = limiter->get_rec();
    
    limiter->calcul(delta_u1+ddu1*0.5,u[NXCELL+1][j]-u[NXCELL][j]);
    du = limiter->get_rec();
    limiter->calcul(delta_v1+ddv1*0.5,v[NXCELL+1][j]-v[NXCELL][j]);
    dv = limiter->get_rec();
    
    h1r[NXCELL][j] = h[NXCELL][j]+dh*0.5;
    h1l[NXCELL][j] = h[NXCELL][j]-dh*0.5;
    
    z1r[NXCELL][j] = z[NXCELL][j]+0.5*(dz_h-dh);
    z1l[NXCELL][j] = z[NXCELL][j]+0.5*(dh-dz_h);
    delzc1[NXCELL][j] = z1r[NXCELL][j]-z1l[NXCELL][j];
    delz1[NXCELL-1][j] = z1l[NXCELL][j]-z1r[NXCELL-1][j];
    delz1[NXCELL][j] = z1l[NXCELL+1][j]-z1r[NXCELL][j];
    
    if (h[NXCELL][j]>0.){
      u1r[NXCELL][j] = u[NXCELL][j]+h1l[NXCELL][j]*du*0.5/h[NXCELL][j];
      u1l[NXCELL][j] = u[NXCELL][j]-h1r[NXCELL][j]*du*0.5/h[NXCELL][j];
      v1r[NXCELL][j] = v[NXCELL][j]+h1l[NXCELL][j]*dv*0.5/h[NXCELL][j];
      v1l[NXCELL][j] = v[NXCELL][j]-h1r[NXCELL][j]*dv*0.5/h[NXCELL][j];
    }
    else{
      u1r[NXCELL][j] = u[NXCELL][j]+du*0.5;
      u1l[NXCELL][j] = u[NXCELL][j]-du*0.5;
      v1r[NXCELL][j] = v[NXCELL][j]+dv*0.5;
      v1l[NXCELL][j] = v[NXCELL][j]-dv*0.5;
    } //end if
    
    h1r[0][j] = h[0][j];
    u1r[0][j] = u[0][j];
    v1r[0][j] = v[0][j];
    h1l[NXCELL+1][j] = h[NXCELL+1][j];
    u1l[NXCELL+1][j] = u[NXCELL+1][j];
    v1l[NXCELL+1][j] = v[NXCELL+1][j];
    
  } //end for
  
  
  for (int i=1 ; i<NXCELL+1 ; i++){
    
    ddh1=0.;
    ddz1=0.;
    ddu1=0.;
    ddv1=0.;
    
    hh1 = h[i][0]-2.*h[i][1]+h[i][2];
    uu1 = u[i][0]-2.*u[i][1]+u[i][2];
    vv1 = v[i][0]-2.*v[i][1]+v[i][2];
    
    delta_h1 = h[i][1]-h[i][0];
    delta_u1 = u[i][1]-u[i][0];
    delta_v1 = v[i][1]-v[i][0];
    
    for (int j=1 ; j<NYCELL ; j++){
      hh2 = h[i][j]-2.*h[i][j+1]+h[i][j+2];
      uu2 = u[i][j]-2.*u[i][j+1]+u[i][j+2];
      vv2 = v[i][j]-2.*v[i][j+1]+v[i][j+2];
      
      limiter->calcul(hh1,hh2);
      ddh2 = AMORTENO*limiter->get_rec();
      limiter->calcul(uu1,uu2);
      ddu2 = AMORTENO*limiter->get_rec();
      limiter->calcul(hh1+som_z2[i][j],hh2+som_z2[i][j+1]);
      ddz2 = AMORTENO*limiter->get_rec();
      limiter->calcul(vv1,vv2);
      ddv2 = AMORTENO*limiter->get_rec();
      
      hh1 = hh2;
      uu1 = uu2;
      vv1 = vv2;
      
      delta_h2 = h[i][j+1]-h[i][j];
      delta_u2 = u[i][j+1]-u[i][j];
      delta_v2 = v[i][j+1]-v[i][j];
      
      limiter->calcul(delta_h1+ddh1*0.5,delta_h2-ddh2*0.5);
      a2 = limiter->get_rec();
      limiter->calcul(delta_h1+delta_z2[i][j-1]+ddz1*0.5,delta_h2+delta_z2[i][j]-ddz2*0.5);
      a4 = limiter->get_rec();
      
      limiter->calcul(delta_h1,delta_h2);
      a1 = limiter->get_rec();
      limiter->calcul(2*MODIFENO*a1,a2);
      dh = limiter->get_rec();
      
      limiter->calcul(delta_h1+delta_z2[i][j-1],delta_h2+delta_z2[i][j]);
      a3 = limiter->get_rec();
      limiter->calcul(2*MODIFENO*a3,a4);
      dz_h = limiter->get_rec();
      
      limiter->calcul(delta_u1+ddu1*0.5,delta_u2-ddu2*0.5);
      du = limiter->get_rec();
      limiter->calcul(delta_v1+ddv1*0.5,delta_v2-ddv2*0.5);
      dv = limiter->get_rec();
      
      delta_h1 = delta_h2;
      delta_u1 = delta_u2;
      delta_v1 = delta_v2;
      
      h2r[i][j] = h[i][j]+dh*0.5;
      h2l[i][j] = h[i][j]-dh*0.5;
      
      z2r[i][j] = z[i][j]+0.5*(dz_h-dh);
      z2l[i][j] = z[i][j]+0.5*(dh-dz_h);
      delzc2[i][j] = z2r[i][j]-z2l[i][j];
      delz2[i][j-1] = z2l[i][j]-z2r[i][j-1];
      
      if (h[i][j]>0.){
        u2r[i][j] = u[i][j]+h2l[i][j]*du*0.5/h[i][j];
        u2l[i][j] = u[i][j]-h2r[i][j]*du*0.5/h[i][j];
        v2r[i][j] = v[i][j]+h2l[i][j]*dv*0.5/h[i][j];
        v2l[i][j] = v[i][j]-h2r[i][j]*dv*0.5/h[i][j];
      }
      else{
        u2r[i][j] = u[i][j]+du*0.5;
        u2l[i][j] = u[i][j]-du*0.5;
        v2r[i][j] = v[i][j]+dv*0.5;
        v2l[i][j] = v[i][j]-dv*0.5;
      } //end if
      
      ddh1=ddh2;
      ddz1=ddz2;
      ddu1=ddu2;
      ddv1=ddv2;
      
    } //end for
    
    limiter->calcul(delta_h1+ddh1*0.5,(h[i][NYCELL+1]-h[i][NYCELL]));
    a2 = limiter->get_rec();
    limiter->calcul(delta_h1+delta_z2[i][NYCELL-1]+ddz1*0.5,h[i][NYCELL+1]-h[i][NYCELL]+delta_z2[i][NYCELL]);
    a4 = limiter->get_rec();
    
    limiter->calcul(delta_h1,h[i][NYCELL+1]-h[i][NYCELL]);
    a1 = limiter->get_rec();
    limiter->calcul(2*MODIFENO*a1,a2);
    dh = limiter->get_rec();
    
    limiter->calcul(delta_h1+delta_z2[i][NYCELL-1],h[i][NYCELL+1]-h[i][NYCELL]+delta_z2[i][NYCELL]);
    a3 = limiter->get_rec();
    limiter->calcul(2*MODIFENO*a3,a4);
    dz_h = limiter->get_rec();
    
    limiter->calcul(delta_u1+ddu1*0.5,u[i][NYCELL+1]-u[i][NYCELL]);
    du = limiter->get_rec();
    limiter->calcul(delta_v1+ddv1*0.5,v[i][NYCELL+1]-v[i][NYCELL]);
    dv = limiter->get_rec();
    
    h2r[i][NYCELL] = h[i][NYCELL]+dh*0.5;
    h2l[i][NYCELL] = h[i][NYCELL]-dh*0.5;
    
    z2r[i][NYCELL] = z[i][NYCELL]+0.5*(dz_h-dh);
    z2l[i][NYCELL] = z[i][NYCELL]+0.5*(dh-dz_h);
    delzc2[i][NYCELL] = z2r[i][NYCELL]-z2l[i][NYCELL];
    delz2[i][NYCELL-1] = z2l[i][NYCELL]-z2r[i][NYCELL-1];
    delz2[i][NYCELL] = z2l[i][NYCELL+1]-z2r[i][NYCELL];
    
    if (h[i][NYCELL]>0.){
      u2r[i][NYCELL] = u[i][NYCELL]+h2l[i][NYCELL]*du*0.5/h[i][NYCELL];
      u2l[i][NYCELL] = u[i][NYCELL]-h2r[i][NYCELL]*du*0.5/h[i][NYCELL];
      v2r[i][NYCELL] = v[i][NYCELL]+h2l[i][NYCELL]*dv*0.5/h[i][NYCELL];
      v2l[i][NYCELL] = v[i][NYCELL]-h2r[i][NYCELL]*dv*0.5/h[i][NYCELL];
    }
    else{
      u2r[i][NYCELL] = u[i][NYCELL]+du*0.5;
      u2l[i][NYCELL] = u[i][NYCELL]-du*0.5;
      v2r[i][NYCELL] = v[i][NYCELL]+dv*0.5;
      v2l[i][NYCELL] = v[i][NYCELL]-dv*0.5;
    } //end if
    
    h2r[i][0] = h[i][0];
    u2r[i][0] = u[i][0];
    v2r[i][0] = v[i][0];
    h2l[i][NYCELL+1] = h[i][NYCELL+1];
    u2l[i][NYCELL+1] = u[i][NYCELL+1];
    v2l[i][NYCELL+1] = v[i][NYCELL+1];
    
  } //end for
  
}

ENO_mod::~ENO_mod(){
  
  if (limiter != NULL){
    delete limiter;
    limiter = NULL;
  }
  
  for (int i=1 ; i<=NXCELL ; i++){
    som_z1[i].clear();
    som_z2[i].clear();
  }
  
  som_z1.clear();
  som_z2.clear();
}

