/**
 * @file f_rusanov.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Rusanov flux
 * @details 
 * Numerical flux: Rusanov formulation.
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

#include "f_rusanov.hpp"

F_Rusanov::F_Rusanov()
{
}

void F_Rusanov::calcul(SCALAR h_L, SCALAR u_L, SCALAR v_L, SCALAR h_R, SCALAR u_R, SCALAR v_R)
{

  /**
   * @details
   * Recall that this is reduced to a one-dimensional computation along the normal of the mesh edge.\n
   * If the water heights on the two sides are small, there is no water.
   * Else, Rusanov formulation is used (see \cite Bouchut04):
   * \f[{\cal F}(U_L,U_R)=\frac{F(U_L)+F(U_R)}{2}-c\frac{U_R-U_L}{2},\f]
   * with \f$ c=\max\left(|u_L|+\sqrt{gh_L},|u_R|+\sqrt{gh_R}\right)\f$,
   * \f$U = {}^t (h,hu,hv)\f$ and \f$ F(U) = {}^t (hu, hu^2+gh^2/2, hv^2) \f$.
   * @param[in] h_L water height at the left of the interface where the flux is calculated.
   * @param[in] u_L velocity (in the x direction) at the left of the interface where the flux is calculated.
   * @param[in] v_L velocity (in the y direction) at the left of the interface where the flux is calculated.
   * @param[in] h_R water height at the right of the interface where the flux is calculated.
   * @param[in] u_R velocity (in the x direction) at the right of the interface where the flux is calculated.
   * @param[in] v_R velocity (in the y direction) at the right of the interface where the flux is calculated.
   * @par Modifies
   * Flux#f1 first componant of the numerical flux.\n
   * Flux#f2 second componant of the numerical flux.\n
   * Flux#f3 third componant of the numerical flux.\n
   * Flux#cfl value of the CFL.
   */

  SCALAR c;
  if (h_L < HE_CA && h_R < HE_CA)
  {
    c = 0.;
    f1 = 0.;
    f2 = 0.;
    f3 = 0.;
    cfl = 0.;
  }
  else
  {
    c = max(fabs(u_L) + sqrt(GRAV * h_L), fabs(u_R) + sqrt(GRAV * h_R));
    SCALAR cd = c * 0.5;
    SCALAR q_R = u_R * h_R;
    SCALAR q_L = u_L * h_L;
    f1 = (q_L + q_R) * 0.5 - cd * (h_R - h_L);
    f2 = ((u_L * q_L) + (GRAV_DEM * h_L * h_L) + (u_R * q_R) + (GRAV_DEM * h_R * h_R)) * 0.5 - cd * (q_R - q_L);
    f3 = (q_L * v_L + q_R * v_R) * 0.5 - cd * (h_R * v_R - h_L * v_L);
    cfl = c; //*tx;
  }
}

F_Rusanov::~F_Rusanov()
{
}
