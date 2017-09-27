/**
 * @file choice_friction.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Choice of friction law
 * @details 
 * From the value of the corresponding parameter,
 * calls the chosen friction law.
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

#include "choice_friction.hpp"

Choice_friction::Choice_friction(Parameters &par)
{

  /**
   * @details
   * Defines the friction law from the value given in the parameters file.
   * @param[in] par parameter, contains all the values from the parameters file.
   */

  switch (par.get_fric())
  {
  case 0:
    fric = new No_Friction(par);
    break;
  case 1:
    fric = new Fr_Manning(par);
    break;
  case 2:
    fric = new Fr_Darcy_Weisbach(par);
    break;
  case 3:
    fric = new Fr_Laminar(par);
    break;
  }
}

void Choice_friction::calcul(const TAB &uold, const TAB &vold, const TAB &hnew, const TAB &q1new, const TAB &q2new, SCALAR dt)
{

  /**
   * @details
   * Calls the calculation of the friction law.
   * @param[in] uold velocity in the first direction at the previous time (\f$n\f$ if you are calculating the \f$n+1\f$th time step).
   * @param[in] vold velocity in the second direction at the previous time (\f$n\f$ if you are calculating the \f$n+1\f$th time step).
   * @param[in] hnew water height after the Shallow-Water computation (without friction).
   * @param[in] q1new discharge in the first direction after the Shallow-Water computation (without friction).
   * @param[in] q2new discharge in the second direction after the Shallow-Water computation (without friction).
   * @param[in] dt time step.
   * @note The friction only affects the discharge.
   */

  fric->calcul(uold, vold, hnew, q1new, q2new, dt);
}

TAB Choice_friction::get_q1mod()
{

  /**
   * @details
   * Calls the function to get the discharge in the first direction modified by the friction term.
   * @return Friction#q1mod discharge in the first direction modified by the friction term.
   */

  return fric->get_q1mod();
}

TAB Choice_friction::get_q2mod()
{

  /**
   * @details
   * Calls the function to get the discharge in the second direction modified by the friction term.
   * @return Friction#q2mod discharge in the second direction modified by the friction term.
   */

  return fric->get_q2mod();
}

void Choice_friction::calculSf(const TAB &h, const TAB &u, const TAB &v)
{

  /**
   * @details
   * Calls the calculation of the explicit friction law.
   * @param[in] h water height.
   * @param[in] u velocity in the first direction.
   * @param[in] v velocity in the second direction.
   * @note This term will be used to compute erosion.
   */

  fric->calculSf(h, u, v);
}

TAB Choice_friction::get_Sf1()
{

  /**
   * @details
   * Calls the function to get the explicit friction term in the first direction.
   * @return Friction#Sf1 explicit friction term in the first direction.
   */

  return fric->get_Sf1();
}

TAB Choice_friction::get_Sf2()
{

  /**
   * @details
   * Calls the function to get the explicit friction term in the second direction.
   * @return Friction#Sf2 explicit friction term in the second direction.
   */

  return fric->get_Sf2();
}

Choice_friction::~Choice_friction()
{
  if (fric != NULL)
  {
    delete fric;
    fric = NULL;
  }
}
