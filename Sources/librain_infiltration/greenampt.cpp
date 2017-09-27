/**
 * @file greenampt.cpp
 * @author Marie Rousseau <M.Rousseau@brgm.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief  Green-Ampt law
 * @details
 * %Infiltration law: bi-layer Green-Ampt.
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

#include "greenampt.hpp"

GreenAmpt::GreenAmpt(Parameters &par) : Infiltration(par)
{

  /**
   * @details Initializes the values for Green-Ampt infiltration.
   * @param[in] par parameter, contains all the values from the parameters file.
   * @warning ***: ERROR: the value at the point ***.
   */

  /*----------------------------------------------------------------------- */
  /* Initializations and allocations. */
  Ks.resize(NXCELL + 1);     // i : 1->NXCELL
  Kc.resize(NXCELL + 1);     // i : 1->NXCELL
  dtheta.resize(NXCELL + 1); // i : 1->NXCELL
  Psi.resize(NXCELL + 1);    // i : 1->NXCELL
  zcrust.resize(NXCELL + 1); // i : 1->NXCELL
  imax.resize(NXCELL + 1);   // i : 1->NXCELL
  for (int i = 0; i < NXCELL + 1; i++)
  {
    Ks[i].resize(NYCELL + 1);     // i : 1->NXCELL
    Kc[i].resize(NYCELL + 1);     // i : 1->NXCELL
    dtheta[i].resize(NYCELL + 1); // i : 1->NXCELL
    Psi[i].resize(NYCELL + 1);    // i : 1->NXCELL
    zcrust[i].resize(NYCELL + 1); // i : 1->NXCELL
    imax[i].resize(NYCELL + 1);   // i : 1->NXCELL
  }

  //**********************************************************************
  //Kc: hydraulic conductivity of the (upper) crust
  //**********************************************************************
  /* Initialization of Kc*/
  if (1 == par.get_Kc_init())
  { //initialization from a file
    par.fill_array(Kc, par.get_KcNameFile());

    /* We verify that Kc is not negative in the case of an initialization from a file */
    for (int i = 1; i < NXCELL + 1; i++)
    {
      for (int j = 1; j < NYCELL + 1; j++)
      {
        if (Kc[i][j] < 0)
        { //if hydraulic conductivity of the crust negative
          cerr << par.get_KcNameFile() << ": ERROR: the value at the point " << (i - 0.5) * DX << "    " << (j - 0.5) * DY << " is " << Kc[i][j] << " ;\n it must be positive." << endl;
          exit(EXIT_FAILURE);
        }
      } //end for j
    }   //end for i
  }
  else
  {
    par.fill_array(Kc, par.get_Kc_coef());
  } //end if

  //**********************************************************************
  //Ks: hydraulic conductivity of the (lower) soil.
  //**********************************************************************
  /* Initialization of Ks*/
  if (1 == par.get_Ks_init())
  { //initialization from a file
    par.fill_array(Ks, par.get_KsNameFile());

    /* We verify that Ks is not negative in the case of an initialization from a file */
    for (int i = 1; i < NXCELL + 1; i++)
    {
      for (int j = 1; j < NYCELL + 1; j++)
      {
        if (Ks[i][j] < 0)
        { //if hydraulic conductivity of the soil negative
          cerr << par.get_KsNameFile() << ": ERROR: the value at the point " << (i - 0.5) * DX << "    " << (j - 0.5) * DY << " is " << Ks[i][j] << " ;\n it must be positive." << endl;
          exit(EXIT_FAILURE);
        }
      } //end for j
    }   //end for i
  }
  else
  { //initialization with a number
    par.fill_array(Ks, par.get_Ks_coef());
  } //end if

  //**********************************************************************
  //dtheta: water content.
  //**********************************************************************
  /* Initialization of dtheta*/
  if (1 == par.get_dtheta_init())
  { //initialization from a file
    par.fill_array(dtheta, par.get_dthetaNameFile());

    /*As dtheta is dimensionless, we verify that its value is between zero and one.*/
    for (int i = 1; i < NXCELL + 1; i++)
    {
      for (int j = 1; j < NYCELL + 1; j++)
      {
        if (dtheta[i][j] < 0)
        { //if water content negative
          cerr << par.get_dthetaNameFile() << ": ERROR: the value at the point " << (i - 0.5) * DX << "    " << (j - 0.5) * DY << " is " << dtheta[i][j] << " ;\n it must be positive." << endl;
          exit(EXIT_FAILURE);
        }
        if (dtheta[i][j] > 1)
        { //if water content are above 1
          cerr << par.get_dthetaNameFile() << ": ERROR: the value at the point " << (i - 0.5) * DX << "    " << (j - 0.5) * DY << " is " << dtheta[i][j] << " ;\n it must be lower than 1." << endl;
          exit(EXIT_FAILURE);
        }
      } //end for j
    }   //end for i
  }
  else
  {
    par.fill_array(dtheta, par.get_dtheta_coef());
  } //end if

  //**********************************************************************
  //Psi: load pressure.
  //**********************************************************************
  /* Initialization of Psi*/
  if (1 == par.get_Psi_init())
  { //initialization from a file
    par.fill_array(Psi, par.get_PsiNameFile());

    /* We verify that Psi is not negative in the case of an initialization from a file */
    for (int i = 1; i < NXCELL + 1; i++)
    {
      for (int j = 1; j < NYCELL + 1; j++)
      {
        if (Psi[i][j] < 0)
        { //if load pressure negative
          cerr << par.get_PsiNameFile() << ": ERROR: the value at the point " << (i - 0.5) * DX << "    " << (j - 0.5) * DY << " is " << Psi[i][j] << " ;\n it must be positive." << endl;
          exit(EXIT_FAILURE);
        }
      } //end for j
    }   //end for i
  }
  else
  { //initialization with a number
    par.fill_array(Psi, par.get_Psi_coef());
  } //end if

  //**********************************************************************
  //zcrust: thickness of the (upper) crust.
  //**********************************************************************
  /* Initialization of zcrust*/
  if (1 == par.get_zcrust_init())
  { //initialization from a file
    par.fill_array(zcrust, par.get_zcrustNameFile());

    /* We verify that zcrust is not negative in the case of an initialization from a file */
    for (int i = 1; i < NXCELL + 1; i++)
    {
      for (int j = 1; j < NYCELL + 1; j++)
      {
        if (zcrust[i][j] < 0)
        { //if thickness of the crust negative
          cerr << par.get_zcrustNameFile() << ": ERROR: the value at the point " << (i - 0.5) * DX << "    " << (j - 0.5) * DY << " is " << zcrust[i][j] << " ;\n it must be positive." << endl;
          exit(EXIT_FAILURE);
        }
      } //end for j
    }   //end for i
  }
  else
  { //initialization with a number
    par.fill_array(zcrust, par.get_zcrust_coef());
  } //end if

  //**********************************************************************
  //imax: Maximun infiltration rate.
  //**********************************************************************
  /* Initialization of imax*/
  if (1 == par.get_imax_init())
  { //initialization from a file.
    par.fill_array(imax, par.get_imaxNameFile());

    /* We verify that imax is not negative in the case of an initialization from a file */
    for (int i = 1; i < NXCELL + 1; i++)
    {
      for (int j = 1; j < NYCELL + 1; j++)
      {
        if (imax[i][j] < 0)
        {
          cerr << par.get_imaxNameFile() << ": ERROR: the value at the point " << (i - 0.5) * DX << "    " << (j - 0.5) * DY << " is " << imax[i][j] << " ;\n it must be positive." << endl;
          exit(EXIT_FAILURE);
        }
      } //end for j
    }   //end for i
  }
  else
  { //initialization with a number
    par.fill_array(imax, par.get_imax_coef());
  } //end if
}

SCALAR GreenAmpt::capacity(const SCALAR h, const SCALAR Vin_tot, const SCALAR dt, const SCALAR Kc, const SCALAR Ks, const SCALAR dtheta, const SCALAR Psi, const SCALAR zcrust)
{

  /**
   * @details
   * the infiltration capacity is given by:
   * \f[ I_{C}= \left\{\begin{array}{lll} 
   *          \displaystyle K_s(1+\frac{Psi+h}{{Z}_{f}}) & \mbox{ if } & zcrust=0 \\[0.5cm]
   *          \displaystyle K_c(1+\frac{Psi+h}{{Z}_{f}}) & \mbox{ if } & Z_{f} \leq {zcrust} \\[0.5cm]
   *          \displaystyle K_{e}(1+\frac{Psi+h}{Z_{f}}) \end{array}\right.,\f]
   * with the effective hydraulic conductivity
   * \f[ K_{e}=\frac{1}{\frac{1}{Ks}*(1-\frac{zcrust*dtheta}{Vin_tot})+zcrust*\frac{dtheta}{Vin_tot}*\frac{1}{Kc}}\f]
   * @param[in] h water height.
   * @param[in] Vin_tot total infiltrated volume.
   * @param[in] dt time step. 
   * @param[in] Kc hydraulic conductivity of the (upper) crust. 
   * @param[in] Ks hydraulic conductivity of the (lower) soil.
   * @param[in] dtheta water content.
   * @param[in] Psi load pressure.
   * @param[in] zcrust thickness of the (upper) crust.
   * @return <em> \b Ic</em>: infiltration capacity.  
   */
  (void)h;  //unused variable
  (void)dt; //unused variable

  if (zcrust < EPSILON)
  {
    // the first layer does not exists, it degenerates into the usual Green-Ampt model
    Ic = Ks + Ks * Psi * dtheta / Vin_tot;
  }
  else
  {
    if (Vin_tot <= zcrust * dtheta + IE_CA) //Zf<=Zc  Ic=Kc(1+(hf+h)/Zf)
    {
      // the effective test is Vin_tot-IE_CA<=zcrust*dtheta
      Ic = Kc + Kc * Psi * dtheta / Vin_tot;
    }
    else
    {
      Ke = 1 / (1 / Ks * (1 - zcrust * dtheta / Vin_tot) + zcrust * dtheta / Vin_tot * 1 / Kc);
      Ic = Ke + Ke * Psi * dtheta / Vin_tot;
    }
  }

  return Ic;
}

void GreenAmpt::calcul(const TAB &h, const TAB &Vin_tot, const SCALAR dt)
{

  /**
   * @details
   * @param[in] h water height.
   * @param[in] Vin_tot total infiltrated volume.
   * @param[in] dt time step.
   * @par Modifies
   * infiltration#hmod modified water height.\n
   * infiltration#Vin total infiltrated volume containing the current time step. 
   */

  for (int i = 1; i < NXCELL + 1; i++)
  {
    for (int j = 1; j < NYCELL + 1; j++)
    {

      /* To avoid infinite values of Ic, we replace Vin_tot by IE_CA
       if Vin_tot is very small (nearly 0).
       This is necessary for the first time step and when there is no
       (or a very small) infiltration. */

      Vin_tot_loc = max(Vin_tot[i][j], IE_CA);

      /* We compute the infiltration capacity, as Vin_tot cannot vanish.
       Note that if (Kc = 0) or (Ks = 0 and zcrust = 0), Ic = 0. */

      Ic = capacity(h[i][j], Vin_tot_loc, dt, Kc[i][j], Ks[i][j], dtheta[i][j], Psi[i][j], zcrust[i][j]);

      /* To avoid an unrealistic infiltration when Ic is large
       due to a small Vin_tot, the maximum value imax can be used. */

      Ic = min(Ic, imax[i][j]);

      /* The infiltrated volume cannot exceed h
       and the water height is modified consequently. */

      Vin[i][j] = min(h[i][j], Ic * dt);
      hmod[i][j] = h[i][j] - Vin[i][j];

      /*The output value is the infiltrated total volume*/
      Vin[i][j] += Vin_tot[i][j];

    } //end for l
  }   //end for i
}

GreenAmpt::~GreenAmpt()
{
}
