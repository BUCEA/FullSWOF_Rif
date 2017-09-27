/**
 * @file rain_read.hpp 
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief  File configuration.
 * @details
 * Initialization of the rain:
 * the values are read in a file.
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

#ifndef RAIN_HPP
#include "rain.hpp"
#endif

#ifndef RAIN_READ_HPP
#define RAIN_READ_HPP

/** @class Rain_read
 * @brief File configuration
 * @details
 * Class that initializes the rain
 * to the values read in a file.
 */

class Rain_read: public Rain{
  
  public:
  
    /** @brief Constructor */
    Rain_read(Parameters &);
    
    /** @brief Performs the initialization */
    void rain_func(SCALAR , TAB &);
  
    /** @brief Destructor */
    virtual ~Rain_read();
  
  private:
    
    string rain_namefile;
    /* Table that contains the intensity of rain (mm/h) */
    vector<SCALAR> rain_tab;
    /* List of the times of changes of intensity */
    vector<SCALAR> rain_times;
    /* Index of the table rain_times */
    int ind;
    int num_lin; //iterator
    SCALAR time_value, intensity_value;
    SCALAR time_value_prev;
    string line;// string to store a line of the input file
    char car;//First character of a commentary line
    
};
#endif
