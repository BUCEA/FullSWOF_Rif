/**
 * @file choice_output.cpp
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Choice of output format
 * @details 
 * From the value of the corresponding parameter,
 * calls the savings in the chosen format.
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

#include "choice_output.hpp"

Choice_output::Choice_output(Parameters & par){
  
  /**
   * @details
   * Defines the output format from the value given in the parameters file.
   * @param[in] par parameter, contains all the values from the parameters file.
   */
  
  switch (par.get_output()){
    case 0:
      out = new No_Evolution_File(par);
      break;
    
    case 1:
      out = new Gnuplot(par);
      break;
    
    case 2:
      out = new Vtk_Out(par);
      break;
  }
}

void Choice_output::write(TAB h, TAB u, TAB v, TAB z, SCALAR time){
  
  /**
   * @details
   * Calls the saving of the current time.
   * @param[in] h water height.
   * @param[in] u first componant of the velocity.
   * @param[in] v second componant of the velocity.
   * @param[in] z topography.
   * @param[in] time value of the current time.
   */
  
  out->write(h,u,v,z,time);
}

void Choice_output:: check_vol(SCALAR time,SCALAR dt, SCALAR Vol_rain_tot,SCALAR  Vol_inf,SCALAR  Vol_of,SCALAR Vol_bound_tot ){
  
  /**
   * @details Calls the saving of the infiltrated and rain volumes.
   * @param[in] time current time.
   * @param[in] dt time step.
   * @param[in] Vol_rain_tot total rain volume.
   * @param[in] Vol_inf volume of infiltrated water.
   * @param[in] Vol_of volume of overland flow.
   * @param[in] Vol_bound_tot total volume of water at the boundary.
   */
  
  out->check_vol(time,dt,Vol_rain_tot, Vol_inf,Vol_of, Vol_bound_tot);

}

void Choice_output:: result(SCALAR time, const clock_t cpu, SCALAR Vol_rain, SCALAR Vol_inf, SCALAR Vol_of, const SCALAR FROUDE, const int NBITER,SCALAR vol_output){
  
  /**
   * @details Calls the saving of the global values.
   * @param[in] time elapsed time.
   * @param[in] cpu CPU time.
   * @param[in] Vol_rain total rain volume.
   * @param[in] Vol_inf total volume of infiltrated water.
   * @param[in] Vol_of total volume of overland flow.
   * @param[in] FROUDE mean Froude number (in space) at the final time.
   * @param[in] NBITER number of time steps.
   * @param[in] vol_output total outflow volume at the boundary.
   */
  
  out->result(time,cpu, Vol_rain, Vol_inf, Vol_of, FROUDE,  NBITER, vol_output);
  
}

void Choice_output::initial(TAB z,TAB h,TAB u,TAB v){
  
  /**
   * @details
   * Calls the saving of the inital time.
   * @param[in] z topography.
   * @param[in] h water height.
   * @param[in] u first componant of the velocity.
   * @param[in] v second componant of the velocity.
   */
  
  out->initial(z, h, u,v);
  
}

void Choice_output::initial_rif(const TAB & rain_c, const TAB & infi_c, const TAB & fric_c) const
{
  /**
   * @details
   * Calls the saving of the initial time.
   * @param[in] rain_c rainfall choice.
   * @param[in] infi_c infiltration choice.
   * @param[in] fric_c friction choice.
  */
  out->initial_rif(rain_c, infi_c, fric_c);
}

void Choice_output::final(TAB z,TAB h,TAB u,TAB v){
  
  /**
   * @details
   * Calls the saving of the final time.
   * @param[in] z topography.
   * @param[in] h water height.
   * @param[in] u first componant of the velocity.
   * @param[in] v second componant of the velocity.
   */
  
  out->final(z, h, u,v);
  
}

SCALAR Choice_output::boundaries_flux(SCALAR time, TAB & flux_u, TAB & flux_v,SCALAR dt,SCALAR dt_first,int ORDER,int verif){
  
  /**
   * @details
   * Calls the saving of the cumulative flux on the boundaries.
   * @param[in] time current time.
   * @param[in] flux_u flux on the left and right boundaries (m^2/s).
   * @param[in] flux_v flux on the bottom and top boundaries (m^2/s).
   * @param[in] dt current time step.
   * @param[in] dt_first previous time step.
   * @param[in] ORDER order of scheme.
   * @param[in] verif parameter to know if we removed the computation with the previous time step (dt_first).
   */
  
  return (out->boundaries_flux( time, flux_u, flux_v, dt, dt_first, ORDER, verif));
}

void Choice_output::boundaries_flux_LR(SCALAR time, TAB LR_flux){
  
  /**
   * @details
   * Calls the saving of the fluxes on the left and right boundaries.
   * @param[in] time current time.
   * @param[in] LR_flux flux on the left and right boundaries (m^2/s).
   */
  
  out->boundaries_flux_LR(time, LR_flux);
}

void Choice_output::boundaries_flux_BT(SCALAR time, TAB BT_flux)
{
  
  /*
   * @details
   * Calls the saving of the fluxes on the top and bottom boundaries.
   * @param[in] time current time.
   * @param[in] BT_flux flux on the bottom and tom boundaries (m^2/s).
   */
  
  out->boundaries_flux_BT(time, BT_flux);
}


Choice_output::~Choice_output(){
  if (out != NULL){
    delete out;
    out = NULL;
  }
}
