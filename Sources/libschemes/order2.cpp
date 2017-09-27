/**
 * @file order2.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Order 2 scheme
 * @details 
 * Numerical scheme: at order 2.
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

#include "order2.hpp"

Order2::Order2(Parameters & par):Scheme(par){
  
  /**
   * @details   
   * Initializations, definition of the reconstruction and creation of 3 vectors for this reconstruction.
   * @param[in] par parameter, contains all the values from the parameters file.
   */
  
  Vin1.resize(NXCELL+1); // i : 1->NXCELL
  Vin2.resize(NXCELL+1); // i : 1->NXCELL
  
  hsa.resize(NXCELL+1); // i : 1->NXCELL
  usa.resize(NXCELL+1); // i : 1->NXCELL
  vsa.resize(NXCELL+1); // i : 1->NXCELL
  qsa1.resize(NXCELL+1); // i : 1->NXCELL
  qsa2.resize(NXCELL+1); // i : 1->NXCELL
  
  for (int i=1 ; i<=NXCELL ; i++){
    Vin1[i].resize(NYCELL+1); // l : 1->NYCELL
    Vin2[i].resize(NYCELL+1); // l : 1->NYCELL
    hsa[i].resize(NYCELL+1); // l : 1->NYCELL
    usa[i].resize(NYCELL+1); // l : 1->NYCELL
    vsa[i].resize(NYCELL+1); // l : 1->NYCELL
    qsa1[i].resize(NYCELL+1); // l : 1->NYCELL
    qsa2[i].resize(NYCELL+1); // l : 1->NYCELL
   
  }
  
  rec = new Choice_reconstruction(par,z);
  
  for (int i=1 ; i<=NXCELL ; i++){
    for (int j=1 ; j<=NYCELL ; j++){
      Vin1[i][j] = 0.;
      Vin2[i][j] = 0.;
    }
  }
}

Order2::~Order2(){
  
  for (int i=1 ; i<=NXCELL ; i++){
    Vin1[i].clear();
    Vin2[i].clear();
    hsa[i].clear();
    usa[i].clear();
    vsa[i].clear();
    qsa1[i].clear();
    qsa2[i].clear();
  }
  
  hsa.clear();
  usa.clear();
  vsa.clear();
  qsa1.clear();
  qsa2.clear();
  
  Vin1.clear();
  Vin2.clear();
  
  if (rec != NULL){
    delete rec;
    rec = NULL;
  }
}

void Order2::calcul(){
  
  /**
   * @details Performs the second order numerical scheme.
   * @note In DEBUG mode, the programme will save another file with volumes of water.
   * @warning order2: WARNING: the computation finished because the maximum number of time steps was reached (see MAX_ITER in misc.hpp)
   */
  
  time(&start);
  //time's iteration beginning 
  
  while (T > cur_time  && n < MAX_ITER+1){
    
    if (1 == verif){
      dt1=dt_max;
      
      // save the data in huz_evolution.dat
      if (cur_time >= T_output){
        out->write(h,u,v,z,cur_time);
        T_output+=dt_output;
      }// end if
      
      //boundary conditions
      boundary(h,u,v,cur_time,NXCELL,NYCELL);
      
      for (int i=1 ; i<NXCELL+1 ; i++){
        for (int j=1 ; j<NYCELL+1 ; j++){
          if (h[i][j]<=HE_CA) {
            h[i][j]=0.;
            u[i][j] = 0.;
            v[i][j] = 0.;
          }
          if (fabs(u[i][j])<=VE_CA){
            u[i][j]=0.;
            q1[i][j]=0.;
          }
          if (fabs(v[i][j])<=VE_CA){
            v[i][j]=0.;
            q2[i][j]=0.;
          }
        } //end for j
      } //end for i
      
      // Reconstruction for order 2
      rec->calcul(h,u,v,z,delzc1,delzc2,delz1,delz2,h1r,u1r,v1r,h1l,u1l,v1l,h2r,u2r,v2r,h2l,u2l,v2l);
      
    }
    else{//if verif==0, ie dt2<dt1
    /*We go back to the values of the previous time */
    for (int i=1 ; i<NXCELL+1 ; i++)
      for (int j=1 ; j<NYCELL+1 ; j++)
        Vin1[i][j]=Vin2[i][j];
    }//end for if verif==1
    
    maincalcflux(CFL_FIX,T, cur_time, dt_max, dt1, dt1);
    
    dt1=min(T-cur_time,dt1);
    
    tx=dt1/DX;
    ty=dt1/DY;
    
    maincalcscheme(h,u,v,q1,q2,hs,us,vs,qs1,qs2,Vin1,cur_time,dt1,n);
    //    dt2=dt_max;
    dt2=dt1;
    
    //boundary conditions
    boundary(hs,us,vs,cur_time+dt1,NXCELL,NYCELL);
    
    for (int i=1 ; i<NXCELL+1 ; i++){
      for (int j=1 ; j<NYCELL+1 ; j++){
        if (hs[i][j]<=HE_CA){
          hs[i][j]=0.;
          us[i][j] = 0.;
          vs[i][j] = 0.;
          qs1[i][j]=0.;
          qs2[i][j]=0.;
        }
        if (fabs(us[i][j])<=VE_CA){
          us[i][j]=0.;
          qs1[i][j]=0.;
        }
        if (fabs(vs[i][j])<=VE_CA){
          vs[i][j]=0.;
          qs2[i][j]=0.;
        }
      } //end for j
    } //end for i
    
    //Reconstruction for order 2 
    rec->calcul(hs,us,vs,z,delzc1,delzc2,delz1,delz2,h1r,u1r,v1r,h1l,u1l,v1l,h2r,u2r,v2r,h2l,u2l,v2l);
    
    //common part
    
    maincalcflux(CFL_FIX,T,cur_time+dt1,dt_max,dt2,dt2);
    
    if (dt2<dt1){
      dt1=dt2;
      tx=dt1/DX;
      ty=dt1/DY;
      verif=0;
    }
    else{
      
      //Added to do calculus at the beginning 
      verif=1;
      maincalcscheme(hs,us,vs,qs1,qs2,hsa,usa,vsa,qsa1,qsa2,Vin1,cur_time,dt1,n+1);
      
      /*the values of height_Vinf_tot and height_of_tot are put to zero
       to compute infiltrated and overland flow volume*/
      height_Vinf_tot=ZERO;
      height_of_tot=ZERO;
      //Heun method (order 2 in time)
      for (int i=1 ; i<NXCELL+1 ; i++){
        for (int j=1 ; j<NYCELL+1 ; j++){
          if (hsa[i][j]<=HE_CA){
            hsa[i][j]=0.;
          }
          if (abs(usa[i][j])<=VE_CA){
            usa[i][j]=0.;
          }
          if (abs(vsa[i][j])<=VE_CA){
            vsa[i][j]=0.;
          }
          tmp = 0.5*(h[i][j]+hsa[i][j]);
          if (tmp>=HE_CA){
            q1[i][j] = 0.5*(h[i][j]*u[i][j]+hsa[i][j]*usa[i][j]);
            u[i][j] = q1[i][j]/tmp;
            q2[i][j] = 0.5*(h[i][j]*v[i][j]+hsa[i][j]*vsa[i][j]);
            v[i][j] = q2[i][j]/tmp;
            h[i][j] = tmp;
          }
          else{
            u[i][j] = 0.;
            q1[i][j] = 0.;
            v[i][j] = 0.;
            q2[i][j] = 0.;
            h[i][j] = 0.;
          }
          Vin_tot[i][j] = (Vin1[i][j] + Vin2[i][j])*0.5;
          
          Vin1[i][j]=Vin_tot[i][j];
          Vin2[i][j]=Vin_tot[i][j];
          
          height_Vinf_tot+= Vin_tot[i][j];
          height_of_tot+=h[i][j];
        } //end for j
      } //end for i
      
      /*Computation of the cumulated infiltrated volume*/
      Vol_inf_tot_cumul=height_Vinf_tot*DX*DY;
      
      /*Computation of the overland flow volume*/
      Vol_of_tot=height_of_tot*DX*DY;
      
      cur_time=cur_time+dt1;
      n=n+1;
      
#ifdef DEBUG
      out->check_vol(cur_time,dt1,Volrain_Tot,Vol_inf_tot_cumul,Vol_of_tot,Total_volume_outflow);
#endif
      
      //Displays the percentage of elapsed time
      cout  << '\r' << '\t' << "[" << int((cur_time/T)*100) << "%] done" ;
      
    }//end for else dt2<dt1
    
  } //end for while : loop in time
  
  //Verifies the reason why the run finished
  if (n>MAX_ITER){
    cerr << "\n order2: WARNING: the computation finished because the maximum number of time steps was reached (see MAX_ITER in misc.hpp)\n";
  }
  
  //Saves the last results
  out->write(h,u,v,z,cur_time);
  
  
  //to give the computing time
  time(&end);
  timecomputation=difftime(end,start);
  cpu_time = clock();
  
  //to inform about the froude number
  Fr=froude_number(hs,us,vs);
  
  
  // The quantity of water outflow is the sum of flux at the boundary multiply by one cell area, so
  //Outflow volum = (fluxNxcell_cum_T+fluxNycell_cum_T)*DX*DY*n1+(fluxx0_cum_T+fluxy0_cum_T)*DX*DY*n2
  //In this case n1=1 and n2=-1 because we consider the direction of the flow
  // is from left to the right (x=0 to x=Nxcell) and from bottom to the top (y=0 to y=NYCELL)  
  out->result(timecomputation,cpu_time,Volrain_Tot, Vol_inf_tot_cumul, Vol_of_tot,Fr,n,Total_volume_outflow);
  //storage of h u v value in the final time
  out->final(z, h, u,v);
  
}
