/**
 * @file choice_flux.hpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.01
 * @date 2016-01-04
 *
 * @brief Choice of numerical flux
 * @details 
 * From the value of the corresponding parameter,
 * calls the chosen numerical flux.
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


#ifndef FLUX_HPP
#include "flux.hpp"
#endif

#ifndef F_RUSANOV_HPP
#include "f_rusanov.hpp"
#endif

#ifndef F_HLL_HPP
#include "f_hll.hpp"
#endif

#ifndef F_HLL2_HPP
#include "f_hll2.hpp"
#endif

#ifndef F_HLLC_HPP
#include "f_hllc.hpp"
#endif

#ifndef F_HLLC2_HPP
#include "f_hllc2.hpp"
#endif

#ifndef CHOICE_FLUX_HPP
#define CHOICE_FLUX_HPP

/** @class Choice_flux
 * @brief Choice of numerical flux
 * @details 
 * Class that calls the numerical flux chosen in the parameters file.
 */


class Choice_flux{
  
  public :
    /** @brief Constructor */
    Choice_flux(int);
  
    /** @brief Calculates the numerical flux */
    void calcul(SCALAR,SCALAR,SCALAR,SCALAR,SCALAR,SCALAR);
  
    /** @brief Sets the variable Flux#tx */
    void set_tx(SCALAR);
  
    /** @brief Gives the first component of the numerical flux */
    SCALAR get_f1();
  
    /** @brief Gives the second component of the numerical flux */
    SCALAR get_f2();
  
    /** @brief Gives the third component of the numerical flux */
    SCALAR get_f3();
  
    /** @brief Gives the CFL value */
    SCALAR get_cfl();
  
    /** @brief Destructor */
    virtual ~ Choice_flux();
    
  private :
    /** The chosen numerical flux. */
    Flux * F;
};
#endif
