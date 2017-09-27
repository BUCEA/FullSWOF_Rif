/**
 * @file flux.hpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Numerical flux
 * @details 
 * Common part for all the numerical fluxes.
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

#include "parameters.hpp"

#ifndef FLUX_HPP
#define FLUX_HPP

/** @class Flux
 * @brief Numerical flux
 * @details 
 * Class that contains all the common declarations for the numerical fluxes.
 */

class Flux{
  
  public :
  
    /** @brief Constructor */
    Flux();
  
    /** @brief Function to be specified in each numerical flux */
    virtual void calcul(SCALAR,SCALAR,SCALAR,SCALAR,SCALAR,SCALAR) =0;
  
    /** @brief Sets the variable Flux#tx */
    void set_tx(SCALAR);
  
    /** @brief Gives the first component of the numerical flux */
    SCALAR get_f1() const;
  
    /** @brief Gives the second component of the numerical flux */
    SCALAR get_f2() const;
  
    /** @brief Gives the third component of the numerical flux */
    SCALAR get_f3() const;
  
    /** @brief Gives the CFL value */
    SCALAR get_cfl() const;
  
    /** @brief Destructor */
    virtual ~Flux();
    
  protected :
    
    /** First component of the numerical flux. */
    SCALAR f1;
    /** Second component of the numerical flux. */
    SCALAR f2;
    /** Third component of the numerical flux. */
    SCALAR f3;
    /** CFL value. */
    SCALAR cfl;
    /** Value of dt/dx. */
    SCALAR tx;
    
};
#endif

