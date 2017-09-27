/**
 * @file huv_generated_radial_dam_dry.hpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Dry radial dam break configuration
 * @details 
 * Initialization of the water height and of the velocity:
 * case of a radial dam break on a dry domain.
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

#ifndef INITIALIZATION_HUV_HPP
#include "initialization_huv.hpp"
#endif

#ifndef HUV_GENERATED_RADIAL_DAM_DRY_HPP
#define HUV_GENERATED_RADIAL_DAM_DRY_HPP

/** @class Huv_generated_Radial_Dam_dry
 * @brief Dry radial dam break configuration
 * @details 
 * Class that initializes the water height and the velocity
 * for a radial dam break on a dry domain.
 */


class Huv_generated_Radial_Dam_dry: public Initialization_huv{
  
  public:
    
    /** @brief Constructor */
    Huv_generated_Radial_Dam_dry(Parameters &);
  
    /** @brief Performs the initialization */
    void initialization(TAB &,TAB &,TAB &);
  
    /** @brief Destructor */
    virtual ~Huv_generated_Radial_Dam_dry();
    
  private:
    
    SCALAR dvalh,dvalhh,dvalu,dvalv;
    int nmid_x,nmid_y,nradius;
    
};
#endif
