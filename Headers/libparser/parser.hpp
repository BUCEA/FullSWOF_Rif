/**
 * @file parser.hpp
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2010-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief %Parser
 * @details 
 * Reads the input file.
 *
 * @copyright License Cecill-V2  \n
 * <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>
 *
 * (c) CNRS - Universite d'Orleans - INRA (France)
 */
/*
 * This file is part of FullSWOF_2D software. 
 * <https://sourcesup.renater.fr/projects/fullswof-1d/> 
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

// When used in FullSWOF_2D:
#include "misc.hpp"

// To be compiled independently:
// comment the previous line and uncomment the 5 following commands
/*
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;
*/

#ifndef PARSER_HPP
#define PARSER_HPP

/** @class Parser
 * @brief %Parser to read the entries
 * @details 
 * Class that reads the input file writen as
 * description \<variable\>:: value # comment
 * and keep the values after the "::" ignoring the comments that begin with a "#".
 */


class Parser {
    
  public:
  
    /** @brief Constructor */
    Parser(const char *);
    
    /** @brief Returns the value of the variable */
    string  getValue(const char *);
    
    /** @brief Destructor */
    virtual ~Parser();
    
  private :
    int i,j, nblines; // nbline = length of data
    string ch;
    string * data; // tabular that contains the strings "description <variable>:: value"
    string value;
    size_t found;
};
#endif
