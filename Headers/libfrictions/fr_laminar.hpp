/**
 * @file fr_laminar.hpp
 * @author Carine Lucas <carine.lucas@univ-orleans.fr> (2014-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief laminar law
 * @details 
 * %Friction law: laminar.
 *
 * @copyright License Cecill-V2 \n
 * <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>
 *
 * (c) CNRS - Universite d'Orleans - INRA (France)
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


#ifndef FRICTION_HPP
#include "friction.hpp"
#endif

#ifndef FR_LAMINAR_HPP
#define FR_LAMINAR_HPP

/** @class Fr_Laminar
 * @brief Laminar law
 * @details 
 * General formulation: \f$ S_f=\nu \frac{1}{gh}\frac{U}{h}\f$.
 * This term is integrated in the code thanks to an implicit method.
 */

class Fr_Laminar: public Friction{
  
  public :
  
    /** @brief Constructor */
    Fr_Laminar(Parameters &);
  
    /** @brief Calculates the laminar friction term */
    void calcul(const TAB &, const TAB &, const TAB &, const TAB &, const TAB &, SCALAR);
    
    /** @brief Calculates the explicit Manning friction term */
    void calculSf(const TAB &, const TAB &, const TAB & );
  
    /** @brief Destructor */
    virtual ~Fr_Laminar();
};
#endif

