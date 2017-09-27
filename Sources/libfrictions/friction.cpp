/**
 * @file friction.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief %Friction law
 * @details 
 * Common part for all the friction laws.
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

#include "friction.hpp"

Friction::Friction(Parameters &par) : NXCELL(par.get_Nxcell()), NYCELL(par.get_Nycell()), DX(par.get_dx()), DY(par.get_dy())
{

  /**
   * @details
   * Defines the number of cells, the space steps and initializes Friction#Fric_tab, Friction#q1mod, Friction#q2mod, Friction#Sf1, Friction#Sf2.
   * @param[in] par parameter, contains all the values from the parameters file.
   * @warning ***: ERROR: the value at the point ***.
   */

  Fric_tab.resize(NXCELL + 1); // i : 1->NXCELL
  for (int i = 1; i < NXCELL + 1; i++)
  {
    Fric_tab[i].resize(NYCELL + 1); // j : 1->NYCELL
  }

  //**********************************************************************
  ///Initialization of Friction
  //**********************************************************************

  if (1 == par.get_fric_init())
  { //initialization from a file
    par.fill_array(Fric_tab, par.get_frictionNameFile());

    /* We verify that the friction is not negative in the case of an initialization from a file */
    for (int i = 1; i < NXCELL + 1; i++)
    {
      for (int j = 1; j < NYCELL + 1; j++)
      {
        if (Fric_tab[i][j] < 0)
        { //if the friction is negative
          cerr << par.get_frictionNameFile() << ": ERROR: the value at the point " << (i - 0.5) * DX << "    " << (j - 0.5) * DY << " is " << Fric_tab[i][j] << " ;\n it must be positive." << endl;
          exit(EXIT_FAILURE);
        }
      } //end for j
    }   //end for i
  }
  else
  {
    par.fill_array(Fric_tab, par.get_friccoef());
  } //end if

  q1mod.resize(NXCELL + 1); // i : 1->NXCELL
  q2mod.resize(NXCELL + 1); // i : 1->NXCELL
  Sf1.resize(NXCELL + 1);   // i : 1->NXCELL
  Sf2.resize(NXCELL + 1);   // i : 1->NXCELL
  for (int i = 1; i < NXCELL + 1; i++)
  {
    q1mod[i].resize(NYCELL + 1); // j : 1->NYCELL
    q2mod[i].resize(NYCELL + 1); // j : 1->NYCELL
    Sf1[i].resize(NYCELL + 1);   // j : 1->NYCELL
    Sf2[i].resize(NYCELL + 1);   // j : 1->NYCELL
  }

  par.fill_array(q1mod, ZERO);
  par.fill_array(q2mod, ZERO);
  par.fill_array(Sf1, ZERO);
  par.fill_array(Sf2, ZERO);
}

TAB Friction::get_q1mod() const
{

  /**
   * @details
   * @return Friction#q1mod discharge in the first direction modified by the friction term.
   */

  return q1mod;
}

TAB Friction::get_q2mod() const
{

  /**
   * @details
   * @return Friction#q2mod discharge in the second direction modified by the friction term.
   */

  return q2mod;
}

TAB Friction::get_Sf1() const
{

  /**
   * @details
   * @return Friction#Sf1 explicit friction term in the first direction.
   */

  return Sf1;
}

TAB Friction::get_Sf2() const
{

  /**
   * @details
   * @return Friction#Sf2 explicit friction term in the second direction.
   */

  return Sf2;
}

Friction::~Friction()
{
  /**
   * @details
   * Deallocation of Friction#Fric_tab, Friction#q1mod, Friction#q2mod, Friction#Sf1, Friction#Sf2
   */

  for (int i = 1; i < NXCELL + 1; i++)
  {
    Fric_tab[i].clear();
    q1mod[i].clear();
    q2mod[i].clear();
    Sf1[i].clear();
    Sf2[i].clear();
  }

  Fric_tab.clear();
  q1mod.clear();
  q2mod.clear();
  Sf1.clear();
  Sf2.clear();
}
