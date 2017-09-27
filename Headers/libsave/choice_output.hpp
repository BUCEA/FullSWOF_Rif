/**
 * @file choice_output.hpp
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief Choice of output format
 * @details 
 * From the value of the corresponding parameter,
 * calls the savings in the chosen format.
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

#ifndef OUTPUT_HPP
#include "output.hpp"
#endif

#ifndef GNUPLOT_HPP
#include "gnuplot.hpp"
#endif

#ifndef VTK_OUT_HPP
#include "vtk_out.hpp"
#endif

#ifndef NO_EVOLUTION_FILE_HPP
#include "no_evolution_file.hpp"
#endif

#ifndef CHOICE_OUTPUT_HPP
#define CHOICE_OUTPUT_HPP

/** @class Choice_output
 * @brief Choice of output format
 * @details 
 * From the value of the corresponding parameter,
 * calls the savings in the chosen format.
 */


class Choice_output{
  
  public :
    
    /** @brief Constructor */
    Choice_output(Parameters &);
    
    /** @brief Save the current time */
    void write(TAB ,TAB ,TAB  ,TAB ,SCALAR );
    
    /** @brief Saves the infiltrated and rain volumes */
    void check_vol(SCALAR,SCALAR,SCALAR,SCALAR ,SCALAR ,SCALAR);
    
    /** @brief Saves global values */
    void result(SCALAR, const clock_t, SCALAR , SCALAR , SCALAR ,  const SCALAR , const int,SCALAR );
    
    /** @brief Saves the initial time */
    void initial(TAB, TAB, TAB,TAB);

    /** @brief saves the initial rainfall, infiltration and friction choice.*/
    void initial_rif(const TAB &, const TAB &, const TAB &) const;
    
    /** @brief Saves the final time */
    void final(TAB, TAB, TAB,TAB);
    
    /** @brief Saves the cumulated fluxes on the boundaries */
    SCALAR boundaries_flux(SCALAR, TAB &, TAB &, SCALAR, SCALAR,int,int);
    
    /** @brief Saves the fluxes on the left and right boundaries */
    void boundaries_flux_LR(SCALAR, TAB);
    
    /** @brief Saves the fluxes on the bottom and top boundaries */
    void boundaries_flux_BT(SCALAR, TAB);
    
    /** @brief Destructor */
    virtual ~Choice_output();
    
  private :
    Output * out;
    
};
#endif
