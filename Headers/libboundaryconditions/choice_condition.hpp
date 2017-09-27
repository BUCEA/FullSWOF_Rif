/**
 * @file choice_condition.hpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Choice of boundary condition
 * @details 
 * From the value of the corresponding parameter,
 * calls the chosen boundary condition.
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


#ifndef BOUNDARY_CONDITION_HPP
#include "boundary_condition.hpp"
#endif

#ifndef BC_IMP_HEIGHT_HPP
#include "bc_imp_height.hpp"
#endif

#ifndef BC_WALL_HPP
#include "bc_wall.hpp"
#endif

#ifndef BC_NEUMANN_HPP
#include "bc_neumann.hpp"
#endif

#ifndef BC_PERIODIC_HPP
#include "bc_periodic.hpp"
#endif

#ifndef BC_IMP_DISCHARGE_HPP
#include "bc_imp_discharge.hpp"
#endif

#ifndef CHOICE_CONDITION_HPP
#define CHOICE_CONDITION_HPP

/** @class Choice_condition
 * @brief Choice of boundary condition
 * @details 
 * Class that calls the boundary condition chosen in the parameters file.
 */



class Choice_condition {
  
  public :
	/** @brief Constructor */
    Choice_condition(int,Parameters & , TAB &,int,int);
  
    /** @brief Calculates the boundary condition */
    void calcul(SCALAR, SCALAR, SCALAR, SCALAR, SCALAR, SCALAR, SCALAR, SCALAR, SCALAR,int,int);
  
    /** @brief Gives the water height on the fictive cell */
    SCALAR get_hbound();
  
    /** @brief Gives the normal velocity of the flow on the fictive cell */
    SCALAR get_unormbound();
  
    /** @brief Gives the tangential velocity of the flow on the fictive cell */
    SCALAR get_utanbound();
  
    /** @brief Destructor */
    virtual ~ Choice_condition();
    
  private :
    
    Boundary_condition * C;
    
};
#endif
