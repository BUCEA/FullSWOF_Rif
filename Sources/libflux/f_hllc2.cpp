/**
 * @file f_hll2.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.01
 * @date 2016-01-04
 *
 * @brief HLL flux
 * @details 
 * Numerical flux: Harten, Lax, van Leer reduced formulation.
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

#include "f_hllc2.hpp"

F_HLLC2::F_HLLC2()
{
}

void F_HLLC2::calcul(SCALAR h_L, SCALAR u_L, SCALAR v_L, SCALAR h_R, SCALAR u_R, SCALAR v_R)
{

  /**
   * @details
   * The HLLC approximate Riemann solver is a modification of the basic HLL scheme to account for the contact and shear waves (see \cite Toro01).\n
   * If the water heights on the two sides are small or \f$c_1 \approx c_2 \approx 0\f$, there is no water.
   * \f[{\cal F}(U_L,U_R)=t_{1}F(U_R)+t_{2}F(U_L)-t_{3}(U_R-U_L),\f]
   * with \f[t_1=\frac{\min(c_2,0)-\min(c_1,0)}{c_2-c_1},\quad t_2=1-t_1,\quad t_3=\frac{c_2|c_1|-c_1|c_2|}{2(c_2-c_1)},\f]
   * \f[c_{1}={\inf\limits_{U=U_L,U_R}}\left({\inf\limits_{j\in\{1,2\}}}|\lambda_{j}(U)|\right) \mbox{ and }
   * c_{2}={\sup\limits_{U=U_L,U_R}}\left({\sup\limits_{j\in\{1,2\}}}|\lambda_{j}(U)|\right),\f]   
   * where \f$\lambda_1(U)=u-\sqrt{gh}\f$ and \f$\lambda_2(U)=u+\sqrt{gh}\f$ are the eigenvalues of the Shallow Water system,
   * \f$U = {}^t (h,hu,hv)\f$ and \f$ F(U) = {}^t (hu, hu^2+gh^2/2, hv^2) \f$.
   * \n If we consider the approximate flux HLL
   * \f[F_{i+\frac{1}{2}}=\left( \begin{array}{lll} F^{1}_{i+\frac{1}{2}}  \\[0.2cm]
   * F^{2}_{i+\frac{1}{2}}\\[0.2cm]
   * F^{3}_{i+\frac{1}{2}} \end{array}\right) \f]

   * \n then to obtain the HLLC solver just add the following expression for the third component
   * \f[F^{3}_{i+\frac{1}{2}}=\left\{\begin{array}{lll} F^{1}_{i+\frac{1}{2}}*V_L & \mbox{if}& 0 \leq u_{*}, \\[0.2cm]
   * F^{1}_{i+\frac{1}{2}}*V_R&\mbox{if}&  u_{*}<0, \end{array}\right.\f]
   * Where \f[\displaystyle u_{*}=\frac{c_{1}h_{R}(u_{R}-c_{2})-c_{2}h_{L}(u_{L}-c_{1})}{h_{R}(u_{R}-c_{2})-h_{L}(u_{L}-c_{1})}\f]

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
   * @note Long double are used locally in the computation to avoid numerical approximations.
   */

  if (h_L < HE_CA && h_R < HE_CA)
  {
    f1 = 0.;
    f2 = 0.;
    f3 = 0.;
    cfl = 0.;
  }
  else
  {
    SCALAR grav_h_L = GRAV * h_L;
    SCALAR grav_h_R = GRAV * h_R;
    SCALAR sqrt_grav_h_L = sqrt(grav_h_L);
    SCALAR sqrt_grav_h_R = sqrt(grav_h_R);
    SCALAR q_R = u_R * h_R;
    SCALAR q_L = u_L * h_L;
    SCALAR c1, c2;

    if (h_L < HE_CA)
    {
      c1 = u_R - 2 * sqrt(GRAV * h_R);
    }
    else
    {
      c1 = min(u_L - sqrt_grav_h_L, u_R - sqrt_grav_h_R); // as u-sqrt(grav_h) <= u+sqrt(grav_h)
    }
    if (h_R < HE_CA)
    {
      c2 = u_L + 2 * sqrt(GRAV * h_L);
    }
    else
    {
      c2 = max(u_L + sqrt_grav_h_L, u_R + sqrt_grav_h_R); // as u+sqrt(grav_h) >= u-sqrt(grav_h)
    }
    SCALAR tmp = 1. / (c2 - c1);
    SCALAR t1 = (min(c2, 0.) - min(c1, 0.)) * tmp;
    SCALAR t2 = 1. - t1;
    SCALAR t3 = (c2 * fabs(c1) - c1 * fabs(c2)) * 0.5 * tmp;
    SCALAR c_star = (c1 * h_R * (u_R - c2) - c2 * h_L * (u_L - c1)) / (h_R * (u_R - c2) - h_L * (u_L - c1));

    f1 = t1 * q_R + t2 * q_L - t3 * (h_R - h_L);
    f2 = t1 * (q_R * u_R + grav_h_R * h_R * 0.5) + t2 * (q_L * u_L + grav_h_L * h_L * 0.5) - t3 * (q_R - q_L);
    if (c_star > EPSILON)
    {
      f3 = f1 * v_L;
    }
    else
    {
      f3 = f1 * v_R;
    }
    cfl = max(fabs(c1), fabs(c2)); //cfl is the velocity to compute the cfl condition max(fabs(c1),fabs(c2))*tx with tx=dt/dx
  }
}

F_HLLC2::~F_HLLC2()
{
}
