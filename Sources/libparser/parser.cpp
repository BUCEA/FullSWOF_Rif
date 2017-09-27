/**
 * @file parser.cpp
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

#include "parser.hpp"

Parser::Parser(const char *FILENAME)
{

  /**
   * @details  
   * Constructor: reads the input parameter and copy the data in a tabular.
   * @param[in] FILENAME name of the paramters file.
   * @warning Impossible to open the *** file.
   * @note If the parameters file cannot be opened, the code will exit with failure termination code.
   */

  i = 0;
  j = 0;
  nblines = 0;

  ifstream entries(FILENAME, ios::in);
  if (!entries)
  {
    cout << "Impossible to open the " << FILENAME << " file\n";
    exit(EXIT_FAILURE);
  }

  // the length of data is fixed equal to the number of lines of the parameter file
  // (even if the comment lines should be eliminated)

  // Find the number of lines of the parameters file

  while (!entries.eof())
  {                       // while the end of the file is not attained
    getline(entries, ch); // read a line
    nblines++;            // nblines is the number of lines of the parameter file
  }

  entries.clear();
  entries.seekg(0, ios::beg); // the file will be read again from the beginning
  //istream& seekg (streampos pos);
  //istream& seekg (streamoff off, ios_base::seekdir way);
  //Sets the position of next character to be extracted from the input stream

  data = new string[nblines];

  // data is a tabular where the useful lines of the input parameter are copied
  // the comments (with a #) and the empty lines are not copied

  while (!entries.eof())
  {                       // while the end of the file is not attained
    getline(entries, ch); // the whole line is copied in the string ch
    found = ch.find("#");
    j = int(found); // j is the position of the first #
    if (ch != "" && j != 0)
    { // if the line is not empty and if it is not a comment line (which begins with a #)
      if (j > 0)
      {                                     // if there is a comment (which begins with a #) after the datas
        ch.erase(ch.begin() + j, ch.end()); // erase the comment
      }
      data[i] = ch;
      i++;
    }
  }

  nblines = i - 1; // number of non-empty lines of data

  entries.close();
}

string Parser::getValue(const char *TAG)
{

  /**
   * @details
   * Return the value corresponding to the tag.
   * @param[in] TAG name of the variable with delimiters. such as the "<Nxcell>"
   * @warning No entry for the variable ***.
   * @warning Bad syntax for ***. The syntax is: description \<variable\>:: value
   * @return Value of the variable as a string
   * @note If the value cannot be read correctly, the code will exit with failure termination code.
   */

  // Find the line of the tabular data that contains the string TAG and copy the line in value

  i = 0;
  j = 0;

  while (j <= 0 && i <= nblines)
  {
    value = data[i];
    found = value.find(TAG);
    j = int(found);
    i++;
  }

  // If the string TAG is not in the parameters file, the program must stop

  if (j <= 0)
  {
    cout << "No entry for the variable" << TAG << endl;
    exit(EXIT_FAILURE);
  }

  // In value, erase the description of the variable

  found = value.find("::"); // j is the position of the ::
  j = int(found);
  if (-1 == j)
  { // if the :: are not found
    cout << "Bad syntax for " << TAG << " . The syntax is: description <variable>:: value " << endl;
    exit(EXIT_FAILURE);
  }

  value.erase(0, j + 2); // erase the description of the variable, before the ::

  value.erase(0, int(value.find_first_not_of(" "))); // erase the white spaces before the string

  value.erase(value.begin() + int(value.find_last_not_of(" ")) + 1, value.end()); // erase the white spaces after the string

  return value;
}

Parser::~Parser()
{
  delete[] data;
}
