/**
 * @file muscl.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief %MUSCL reconstruction
 * @details 
 * Linear reconstruction: %MUSCL.
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

#include "muscl.hpp"

MUSCL::MUSCL(Parameters &par, TAB &z) : Reconstruction(par, z)
{

  /**
   * @details
   * Initializations.
   * @param[in] par parameter, contains all the values from the parameters file.
   * @param[in] z topography.
   */
}

void MUSCL::calcul(TAB &h, TAB &u, TAB &v, TAB &z, TAB &delzc1, TAB &delzc2, TAB &delz1, TAB &delz2, TAB &h1r, TAB &u1r, TAB &v1r, TAB &h1l, TAB &u1l, TAB &v1l, TAB &h2r, TAB &u2r, TAB &v2r, TAB &h2l, TAB &u2l, TAB &v2l)
{

  /**
   * @details
   * Calls the calculation of the second order reconstruction in space with %MUSCL formulation, see \cite vanLeer79 \cite Bouchut07b.
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

  for (int j = 1; j < NYCELL + 1; j++)
  {

    delta_h1 = h[1][j] - h[0][j];
    delta_u1 = u[1][j] - u[0][j];
    delta_v1 = v[1][j] - v[0][j];

    for (int i = 1; i < NXCELL + 1; i++)
    {

      delta_h2 = h[i + 1][j] - h[i][j];
      delta_u2 = u[i + 1][j] - u[i][j];
      delta_v2 = v[i + 1][j] - v[i][j];

      limiter->calcul(delta_h1, delta_h2);
      a2 = limiter->get_rec();
      limiter->calcul(delta_h1 + delta_z1[i - 1][j], delta_h2 + delta_z1[i][j]);
      a4 = limiter->get_rec();
      dh = a2;
      dz_h = a4;

      limiter->calcul(delta_u1, delta_u2);
      du = limiter->get_rec();
      limiter->calcul(delta_v1, delta_v2);
      dv = limiter->get_rec();
      h1r[i][j] = h[i][j] + dh * 0.5;
      h1l[i][j] = h[i][j] - dh * 0.5;

      z1r[i][j] = z[i][j] + 0.5 * (dz_h - dh);
      z1l[i][j] = z[i][j] + 0.5 * (dh - dz_h);
      // Long double are used locally in the computation to avoid numerical approximations
      delzc1[i][j] = (SCALAR)((long double)z1r[i][j] - (long double)z1l[i][j]);
      delz1[i - 1][j] = z1l[i][j] - z1r[i - 1][j];

      if (h[i][j] > 0.)
      {
        u1r[i][j] = u[i][j] + h1l[i][j] * du * 0.5 / h[i][j];
        u1l[i][j] = u[i][j] - h1r[i][j] * du * 0.5 / h[i][j];
        v1r[i][j] = v[i][j] + h1l[i][j] * dv * 0.5 / h[i][j];
        v1l[i][j] = v[i][j] - h1r[i][j] * dv * 0.5 / h[i][j];
      }
      else
      {
        u1r[i][j] = u[i][j] + du * 0.5;
        u1l[i][j] = u[i][j] - du * 0.5;
        v1r[i][j] = v[i][j] + dv * 0.5;
        v1l[i][j] = v[i][j] - dv * 0.5;
      } //end if

      delta_h1 = delta_h2;
      delta_u1 = delta_u2;
      delta_v1 = delta_v2;

    } //end for i

    h1r[0][j] = h[0][j];
    u1r[0][j] = u[0][j];
    v1r[0][j] = v[0][j];
    h1l[NXCELL + 1][j] = h[NXCELL + 1][j];
    u1l[NXCELL + 1][j] = u[NXCELL + 1][j];
    v1l[NXCELL + 1][j] = v[NXCELL + 1][j];

    delz1[NXCELL][j] = z1l[NXCELL + 1][j] - z1r[NXCELL][j];
  } //end for j

  for (int i = 1; i < NXCELL + 1; i++)
  {

    delta_h1 = h[i][1] - h[i][0];
    delta_u1 = u[i][1] - u[i][0];
    delta_v1 = v[i][1] - v[i][0];

    for (int j = 1; j < NYCELL + 1; j++)
    {

      delta_h2 = h[i][j + 1] - h[i][j];
      delta_u2 = u[i][j + 1] - u[i][j];
      delta_v2 = v[i][j + 1] - v[i][j];

      limiter->calcul(delta_h1, delta_h2);
      a2 = limiter->get_rec();
      limiter->calcul(delta_h1 + delta_z2[i][j - 1], delta_h2 + delta_z2[i][j]);
      a4 = limiter->get_rec();

      dh = a2;
      dz_h = a4;

      limiter->calcul(delta_u1, delta_u2);
      du = limiter->get_rec();
      limiter->calcul(delta_v1, delta_v2);
      dv = limiter->get_rec();

      h2r[i][j] = h[i][j] + dh * 0.5;
      h2l[i][j] = h[i][j] - dh * 0.5;

      z2r[i][j] = z[i][j] + 0.5 * (dz_h - dh);
      z2l[i][j] = z[i][j] + 0.5 * (dh - dz_h);
      // Long double are used locally in the computation to avoid numerical approximations
      delzc2[i][j] = (SCALAR)((long double)z2r[i][j] - (long double)z2l[i][j]);
      delz2[i][j - 1] = z2l[i][j] - z2r[i][j - 1];

      if (h[i][j] > 0.)
      {
        u2r[i][j] = u[i][j] + h2l[i][j] * du * 0.5 / h[i][j];
        u2l[i][j] = u[i][j] - h2r[i][j] * du * 0.5 / h[i][j];
        v2r[i][j] = v[i][j] + h2l[i][j] * dv * 0.5 / h[i][j];
        v2l[i][j] = v[i][j] - h2r[i][j] * dv * 0.5 / h[i][j];
      }
      else
      {
        u2r[i][j] = u[i][j] + du * 0.5;
        u2l[i][j] = u[i][j] - du * 0.5;
        v2r[i][j] = v[i][j] + dv * 0.5;
        v2l[i][j] = v[i][j] - dv * 0.5;
      } //end if

      delta_h1 = delta_h2;
      delta_u1 = delta_u2;
      delta_v1 = delta_v2;

    } //end for j

    h2r[i][0] = h[i][0];
    u2r[i][0] = u[i][0];
    v2r[i][0] = v[i][0];
    h2l[i][NYCELL + 1] = h[i][NYCELL + 1];
    u2l[i][NYCELL + 1] = u[i][NYCELL + 1];
    v2l[i][NYCELL + 1] = v[i][NYCELL + 1];

    delz2[i][NYCELL] = z2l[i][NYCELL + 1] - z2r[i][NYCELL];
  } //end for i
}

MUSCL::~MUSCL()
{
}
