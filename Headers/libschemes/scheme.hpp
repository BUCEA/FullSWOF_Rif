/**
 * @file scheme.hpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Numerical scheme
 * @details
 * Common part for all the numerical schemes.
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

#ifndef PARAMETERS_HPP
#include "parameters.hpp"
#endif

#ifndef HYDROSTATIC_HPP
#include "hydrostatic.hpp"
#endif

#ifndef CHOICE_CONDITION_HPP
#include "choice_condition.hpp"
#endif

#ifndef CHOICE_FLUX_HPP
#include "choice_flux.hpp"
#endif

#ifndef CHOICE_FRICTION_HPP
#include "choice_friction.hpp"
#endif

#ifndef CHOICE_INFILTRATION_HPP
#include "choice_infiltration.hpp"
#endif

#ifndef CHOICE_INIT_TOPO_HPP
#include "choice_init_topo.hpp"
#endif

#ifndef CHOICE_INIT_HUV_HPP
#include "choice_init_huv.hpp"
#endif

#ifndef CHOICE_INIT_RIF_HPP
#include "choice_init_rif.hpp"
#endif

#ifndef CHOICE_RAIN_HPP
#include "choice_rain.hpp"
#endif

#ifndef CHOICE_OUTPUT_HPP
#include "choice_output.hpp"
#endif

#ifndef CHOICE_RECONSTRUCTION
#include "choice_reconstruction.hpp"
#endif

#ifndef SCHEME_HPP
#define SCHEME_HPP

/** @class Scheme
 * @brief Numerical scheme
 * @details 
 * Class that contains all the common declarations for the numerical schemes.
 */

class Scheme
{

public:
  /** @brief Constructor */
  Scheme(Parameters &);

  /** @brief Function to be specified in each numerical scheme */
  virtual void calcul() = 0;

  /** @brief Allocation of spatialized variables */
  void allocation();

  /** @brief Deallocation of variables */
  void deallocation();

  /** @brief Main calculation of the flux */
  void maincalcflux(SCALAR, SCALAR, SCALAR, SCALAR, SCALAR, SCALAR &);

  /** @brief Main calculation of the scheme */
  void maincalcscheme(TAB &, TAB &, TAB &, TAB &, TAB &, TAB &, TAB &, TAB &, TAB &, TAB &, TAB &, SCALAR, SCALAR, int);

  /** @brief Calls the boundary conditions and affects the boundary values */
  void boundary(TAB &, TAB &, TAB &, SCALAR, const int, const int);

  /** @brief Returns the Froude number */
  SCALAR froude_number(TAB, TAB, TAB);

  /** @brief Destructor*/
  virtual ~Scheme();

protected:
  /** Number of cells in space in the first (x) direction. */
  const int NXCELL;
  /** Number of cells in space in the second (y) direction. */
  const int NYCELL;
  /** Order of the numerical scheme. */
  const int ORDER;
  /** Final time.*/
  const SCALAR T;
  /** Number of times saved. */
  const int NBTIMES;
  /** Type of scheme (fixed cfl or time step). */
  const int SCHEME_TYPE;
  /** Space step in the first (x) direction.*/
  const SCALAR DX;
  /** Space step in the second (y) direction.*/
  const SCALAR DY;
  /** Value of the fixed cfl.*/
  const SCALAR CFL_FIX;
  /** Value of the fixed time step.*/
  SCALAR DT_FIX;
  /** Ratio dt/dx.*/
  SCALAR tx;
  /** Ratio dt/dy.*/
  SCALAR ty;
  /** Time to save the data (evolution file).*/
  SCALAR T_output;
  /** Time step to save the data (evolution file).*/
  SCALAR dt_output;
  /** Friction coefficient. */
  const SCALAR FRICCOEF;
  /** Imposed discharge on the left boundary.*/
  const SCALAR L_IMP_Q;
  /** Imposed water height on the left boundary.*/
  const SCALAR L_IMP_H;
  /** Imposed discharge on the right boundary.*/
  const SCALAR R_IMP_Q;
  /** Imposed water height on the right boundary.*/
  const SCALAR R_IMP_H;
  /** Imposed discharge on the bottom boundary.*/
  const SCALAR B_IMP_Q;
  /** Imposed water height on the bottom boundary.*/
  const SCALAR B_IMP_H;
  /** Imposed discharge on the top boundary.*/
  const SCALAR T_IMP_Q;
  /** Imposed water height on the top boundary.*/
  const SCALAR T_IMP_H;
  /** Topography.*/
  TAB z;
  /** Water height.*/
  TAB h;
  /** First componant of the velocity.*/
  TAB u;
  /** Second componant of the velocity.*/
  TAB v;
  /** Rainfall choice of each cell.*/
  TAB rain_c;
  /** Infiltration choice of each cell.*/
  TAB infi_c;
  /** Friction choice of each cell.*/
  TAB fric_c;
  /** First componant of the discharge.*/
  TAB q1;
  /** Second componant of the discharge.*/
  TAB q2;
  /** Water height after one step of the scheme.*/
  TAB hs;
  /** First componant of the velocity after one step of the scheme.*/
  TAB us;
  /** Second componant of the velocity after one step of the scheme.*/
  TAB vs;
  /** First componant of the discharge after one step of the scheme.*/
  TAB qs1;
  /** Second componant of the discharge after one step of the scheme.*/
  TAB qs2;
  /** First component of the numerical flux along x. */
  TAB f1;
  /** Second component of the numerical flux along x. */
  TAB f2;
  /** Third component of the numerical flux along x. */
  TAB f3;
  /** First component of the numerical flux along y. */
  TAB g1;
  /** Second component of the numerical flux along y. */
  TAB g2;
  /** Third component of the numerical flux along y. */
  TAB g3;
  /** %Hydrostatic reconstruction on the left along x.*/
  TAB h1left;
  /** %Hydrostatic reconstruction on the right along x.*/
  TAB h1right;
  /** %Hydrostatic reconstruction on the left along y.*/
  TAB h2left;
  /** %Hydrostatic reconstruction on the right along y.*/
  TAB h2right;
  /** Variations of the topography along x.*/
  TAB delz1;
  /** Variations of the topography along y.*/
  TAB delz2;
  /** Difference between the reconstructed topographies on the left and on the right boundary of a cell along x.*/
  TAB delzc1;
  /** Difference between the reconstructed topographies on the left and on the right boundary of a cell along y.*/
  TAB delzc2;
  /** Water height on the cell located at the right of the boundary along x.*/
  TAB h1r;
  /** First componant of the velocity on the cell located at the right of the boundary along x.*/
  TAB u1r;
  /** Second componant of the velocity on the cell located at the right of the boundary along x.*/
  TAB v1r;
  /** Water height on the cell located at the left of the boundary along x.*/
  TAB h1l;
  /** First componant of the velocity on the cell located at the left of the boundary along x.*/
  TAB u1l;
  /** Second componant of the velocity on the cell located at the left of the boundary along x.*/
  TAB v1l;
  /** Water height on the cell located at the right of the boundary along y.*/
  TAB h2r;
  /** First componant of the velocity on the cell located at the right of the boundary along y.*/
  TAB u2r;
  /** Second componant of the velocity on the cell located at the right of the boundary along y.*/
  TAB v2r;
  /** Water height on the cell located at the left of the boundary along y.*/
  TAB h2l;
  /** First componant of the velocity on the cell located at the left of the boundary along y.*/
  TAB u2l;
  /** Second componant of the velocity on the cell located at the left of the boundary along y.*/
  TAB v2l;
  /** %Rain intensity.*/
  TAB Rain;
  /** Cumulative volume of infiltrated water (at each point).*/
  TAB Vin_tot;
  /** Beginning of timer.*/
  time_t start;
  /** End of timer.*/
  time_t end;
  /** Duration of the computation.*/
  SCALAR timecomputation;
  /** CPU time.*/
  clock_t cpu_time;
  /** Iterator for the loop in time.*/
  int n;
  /** Mean Froude number.*/
  SCALAR Fr;
  /** Time step in case of fixed cfl. */
  SCALAR dt1;
  /** Maximum time step in case of fixed cfl. */
  SCALAR dt_max;
  /** The current simulation time. */
  SCALAR cur_time;
  /** Space step in the first step in the method Heun.*/
  SCALAR dt_first;
  /** Cumulative volume of rain on the whole domain.*/
  SCALAR Volrain_Tot;
  /** Cumulative outflow volume at the boundary. */
  SCALAR Total_volume_outflow;
  /** Cumulative water height on the whole domain*/
  SCALAR height_of_tot;
  /** Cumulative height of infiltrated water on the whole domain*/
  SCALAR height_Vinf_tot;
  /** Cumulative volume of water infiltrated.*/
  SCALAR Vol_inf_tot_cumul;
  /** Cumulative streammed volume.*/
  SCALAR Vol_of_tot;
  /** The choice of the left boundary condition.*/
  Choice_condition *Lbound;
  /** The choice of the right boundary condition.*/
  Choice_condition *Rbound;
  /** The choice of the bottom boundary condition.*/
  Choice_condition *Bbound;
  /** The choice of the top boundary condition.*/
  Choice_condition *Tbound;
  /** The choice of output.*/
  Choice_output *out;
  /** Flag for the time step.*/
  int verif;

private:
  /** The choice of the left boundary condition.*/
  int nchoice_Lbound;
  /** The choice of the right boundary condition.*/
  int nchoice_Rbound;
  /** The choice of the bottom boundary condition.*/
  int nchoice_Bbound;
  /** The choice of the top boundary condition.*/
  int nchoice_Tbound;
  /** The choice of rain.*/
  Choice_rain *Prain;
  Hydrostatic rec_hydro;
  Choice_flux *flux_num;
  Choice_friction *fric;
  Choice_infiltration *I;
  Choice_init_topo *topo;
  Choice_init_huv *huv_init;
  Choice_init_rif *rif_init;
};
#endif
