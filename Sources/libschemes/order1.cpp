/**
 * @file order1.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Order 1 scheme
 * @details 
 * Numerical scheme: at order 1.
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

#include "order1.hpp"

Order1::Order1(Parameters &par) : Scheme(par)
{

  /**
   * @details
   * @param[in] par parameter, contains all the values from the parameters file.
   */
}

void Order1::calcul()
{

  /**
   * @details Performs the first order numerical scheme.
   * @note In DEBUG mode, the programme will save four other files with boundary fluxes and volumes of water.
   * @warning order1: WARNING: the computation finished because the maximum number of time steps was reached (see MAX_ITER in misc.hpp)
   */

  //boundary conditions
  boundary(h, u, v, cur_time, NXCELL, NYCELL);

  //initialization delz and delzc
  for (int i = 1; i <= NXCELL + 1; i++)
  {
    for (int j = 1; j <= NYCELL; j++)
    {
      delz1[i - 1][j] = z[i][j] - z[i - 1][j];  //Variations of the topography along x
    } //end for j
  }   //end for i

  for (int i = 1; i <= NXCELL; i++)
  {
    for (int j = 1; j <= NYCELL + 1; j++)
    {
      delz2[i][j - 1] = z[i][j] - z[i][j - 1];  //Variations of the topography along y
    } //end for j
  }   //end for i

  for (int i = 1; i <= NXCELL; i++)
  {
    for (int j = 1; j <= NYCELL; j++)
    {
      delzc2[i][j] = 0.;  //Difference between the reconstructed topography on the left and right boundary along y
      delzc1[i][j] = 0.;  //Difference between the reconstructed topography on the left and right boundary along x
    }
  }

  time(&start);

  //time iteration's beginning
  while (T > cur_time && n < MAX_ITER + 1)
  {

    // save the data in huz_evolution.dat
    if (cur_time >= T_output)
    {
      out->write(h, u, v, z, cur_time);
      T_output += dt_output;
    } // end if

    dt1 = dt_max;

    //boundary conditions
    boundary(h, u, v, cur_time, NXCELL, NYCELL);

    for (int i = 1; i < NXCELL + 1; i++)
    {
      for (int j = 1; j < NYCELL + 1; j++)
      {
        if (h[i][j] <= HE_CA)
        {
          h[i][j] = 0.;
          u[i][j] = 0.;
          v[i][j] = 0.;
          q1[i][j] = 0.;
          q2[i][j] = 0.;
        }
        if (fabs(u[i][j]) <= VE_CA)
        {
          u[i][j] = 0.;
          q1[i][j] = 0.;
        }
        if (fabs(v[i][j]) <= VE_CA)
        {
          v[i][j] = 0.;
          q2[i][j] = 0.;
        }
      } //end for j
    }   //end for i

    for (int i = 0; i < NXCELL + 1; i++)
    {
      for (int j = 1; j < NYCELL + 1; j++)
      {
        h1r[i][j] = h[i][j];
        u1r[i][j] = u[i][j];
        v1r[i][j] = v[i][j];
        h1l[i + 1][j] = h[i + 1][j];
        u1l[i + 1][j] = u[i + 1][j];
        v1l[i + 1][j] = v[i + 1][j];
      } //end for j
    }   //end for i

    for (int i = 1; i < NXCELL + 1; i++)
    {
      for (int j = 0; j < NYCELL + 1; j++)
      {
        h2r[i][j] = h[i][j];
        u2r[i][j] = u[i][j];
        v2r[i][j] = v[i][j];
        h2l[i][j + 1] = h[i][j + 1];
        u2l[i][j + 1] = u[i][j + 1];
        v2l[i][j + 1] = v[i][j + 1];
      } //end for j
    }   //end for i

    maincalcflux(CFL_FIX, T, cur_time, dt_max, dt1, dt1);

    dt1 = min(T - cur_time, dt1);

    tx = dt1 / DX;
    ty = dt1 / DY;

    maincalcscheme(h, u, v, q1, q2, hs, us, vs, qs1, qs2, Vin_tot, cur_time, dt1, n);

    for (int i = 1; i < NXCELL + 1; i++)
    {
      for (int j = 1; j < NYCELL + 1; j++)
      {
        h[i][j] = hs[i][j];
        u[i][j] = us[i][j];
        v[i][j] = vs[i][j];
        q1[i][j] = h[i][j] * u[i][j];
        q2[i][j] = h[i][j] * v[i][j];
      } //end for j
    }   //end for i

    /*the values of height_Vinf_tot and height_of_tot are put to zero
     to compute infiltrated and overland flow volume*/
    height_Vinf_tot = ZERO;
    height_of_tot = ZERO;

    for (int i = 1; i < NXCELL + 1; i++)
    {
      for (int j = 1; j < NYCELL + 1; j++)
      {
        if (h[i][j] <= HE_CA)
        {
          h[i][j] = 0.;
          u[i][j] = 0.;
          v[i][j] = 0.;
          q1[i][j] = 0.;
          q2[i][j] = 0.;
        }
        if (fabs(u[i][j]) <= VE_CA)
        {
          u[i][j] = 0.;
          q1[i][j] = 0.;
        }
        if (fabs(v[i][j]) <= VE_CA)
        {
          v[i][j] = 0.;
          q2[i][j] = 0.;
        }

        height_Vinf_tot += Vin_tot[i][j];
        height_of_tot += h[i][j];

      } //end for j
    }   //end for i

    Vol_inf_tot_cumul = height_Vinf_tot * DX * DY;

    Vol_of_tot = height_of_tot * DX * DY;

    cur_time = cur_time + dt1;
    n = n + 1;

#ifdef DEBUG
    out->check_vol(cur_time, dt1, Volrain_Tot, Vol_inf_tot_cumul, Vol_of_tot, Total_volume_outflow);
#endif

    //Displays the percentage of elapsed time
    cout << '\r' << '\t' << "[" << int((cur_time / T) * 100) << "%] done";

  } //end for n : loop in time

  //Verifies the reason why the run finished
  if (n > MAX_ITER)
  {
    cerr << "\n order1: WARNING: the computation finished because the maximum number of time steps was reached (see MAX_ITER in misc.hpp)\n";
  }

  //Saves the last results
  out->write(h, u, v, z, cur_time);

  //to give the computing time
  time(&end);
  timecomputation = difftime(end, start);
  cpu_time = clock();

  //to inform about the froude number
  Fr = froude_number(hs, us, vs);

  // The quantity of water outflow is the sum of flux at the boundary multiply by one cell area, so
  //Outflow volum = (fluxy0_cum_T+fluxNycell_cum_T+fluxx0_cum_T+fluxNxcell_cum_T)*DX*DY
  //In this case n1=1 and n2=-1 because we consider the direction of the flow
  // is from left to the right (x=0 to x=Nxcell) and from bottom to the top (y=0 to y=NYCELL)
  out->result(timecomputation, cpu_time, Volrain_Tot, Vol_inf_tot_cumul, Vol_of_tot, Fr, n, Total_volume_outflow);

  //storage of h u v value in the final time
  out->final(z, h, u, v);
}

Order1::~Order1()
{
}
