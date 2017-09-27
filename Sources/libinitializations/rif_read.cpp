/**
 * @file rif_read.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.00
 * @date 2015-02-19
 *
 * @brief File configuration
 * @details 
 * Initialization of the rainfall, infiltration and friction choice:
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

#include "rif_read.hpp"

Rif_read::Rif_read(Parameters &par) : Initialization_rif(par)
{

  /**
   * @details   
   * Defines the name of the file for the initialization.
   * @param[in] par parameter, contains all the values from the parameters file.
   */

  rif_namefile = par.get_rifNameFile();
}

void Rif_read::initialization(TAB &rain_c, TAB &infi_c, TAB &fric_c)
{

  /**
   * @details
   * Initializes the rainfall, infiltration and friction choice to the values read in the corresponding file.
   * @param[in, out] rain_c rainfall choice.
   * @param[in, out] infi_c infiltration choice.
   * @param[in, out] fric_c friction choice.
   * @warning (rif_namefile): ERROR: cannot open the rif file.
   * @warning (rif_namefile): ERROR: the number of data in this file is too big
   * @warning (rif_namefile): ERROR: line ***.
   * @warning (rif_namefile): WARNING: line ***.
   * @warning (rif_namefile): ERROR: the number of data in this file is too small
   * @warning (rif_namefile): ERROR: the value for the point x *** y *** is missing
   * @note If the file cannot be opened or if the data are not correct, the code will exit with failure termination code.
   */

  SCALAR x, y, r_init, i_init, f_init;
  int row, column;
  int it = 0;      //iterator
  int num_lin = 0; //iterator
  string line;     // string to store a line of the input file
  char car;        //First character of a commentary line

  //Opening the file and verification if it exists.
  ifstream getdata(rif_namefile.c_str(), ios::in);
  if (!getdata)
  {
    //The input file cannot be opened.
    cerr << rif_namefile << ": ERROR: cannot open the rif file\n";
    exit(EXIT_FAILURE);
  }

  //Initialization of rainfall choice to the largest finite representable floating-point number.
  //So that we will be able to check that the user fills the variables r, i and f properly.
  for (int i=1 ; i<NXCELL+1 ; i++){
    for (int j=1 ; j<NYCELL+1 ; j++){
      rain_c[i][j]=MAX_SCAL;
    } //end for j
  } //end for i

  // As long as the end of the file is not reached, the next line is read.
  while (!getdata.eof())
  {
    num_lin++;
    getline(getdata, line, '\n'); // read a line
    istringstream entree(line);
    if (entree >> x >> y >> r_init >> i_init >> f_init)
    {
      //Check that the input file contains the expected number of data.
      if (++it > NXCELL * NYCELL)
      {
        cerr << rif_namefile << ": ERROR: the number of data in this file is too big!" << endl;
        exit(EXIT_FAILURE);
      }
      //We compute the index of the r, i and f arrays from the space variable in the input file.
      row = (int)ceil(x / DX);
      column = (int)ceil(y / DY); //The index corresponds to the smallest integer superior or equal to x/DX or y/DY.

      //Error if the x or y value of the input file is out of the domain.
      if (x < 0)
      {
        cerr << rif_namefile << ": ERROR: line " << num_lin << ", x = " << x << " must be positive." << endl;
        exit(EXIT_FAILURE);
      }
      if (y < 0)
      {
        cerr << rif_namefile << ": ERROR: line " << num_lin << ", y = " << y << " must be positive." << endl;
        exit(EXIT_FAILURE);
      }
      if (x > NXCELL * DX)
      {
        cerr << rif_namefile << ": ERROR: line " << num_lin << ", x = " << x << " must be lower than " << NXCELL * DX << "." << endl;
        exit(EXIT_FAILURE);
      }
      if (y > NYCELL * DY)
      {
        cerr << rif_namefile << ": ERROR: at line " << num_lin << ", y = " << y << " must be lower than " << NYCELL * DY << "." << endl;
        exit(EXIT_FAILURE);
      }

      //Error if the x or y value of the input file is not near the center of the cell.
      if (fabs(x - (row - 0.5) * DX) > DX * RATIO_CLOSE_CELL)
      {
        cout << rif_namefile << ": ERROR: line " << num_lin << "; x = " << x << ";\n This value is not close enough to the center of the cell.\n You may want to replace it with " << (row - 0.5) * DX << endl;
        exit(EXIT_FAILURE);
      }
      if (fabs(y - (column - 0.5) * DY) > DY * RATIO_CLOSE_CELL)
      {
        cout << "column =" << column << endl;
        cout << rif_namefile << ": ERROR: line " << num_lin << "; y = " << y << ";\n This value is not close enough to the center of the cell.\n You may want to replace it with " << (column - 0.5) * DY << endl;
        exit(EXIT_FAILURE);
      }

      // Error if the rainfall choice r_init is negative.
      if (r_init < 0)
      {
        cerr << rif_namefile << ": ERROR: line " << num_lin << "; A NEGATIVE initial rainfall choice was encountered in this file." << endl;
        exit(EXIT_FAILURE);
      }

      //Store the input values into the rain, infi and fric arrays.

      rain_c[row][column] = r_init;
      infi_c[row][column] = i_init;
      fric_c[row][column] = f_init;
    }
    else
    {
      car = '#'; //Initialization of the character used to identify the beginning of a comment
      istringstream entree_car(line);
      entree_car >> car;
      if (car != '#')
      {
        cout << rif_namefile << ": WARNING: line " << num_lin << "; a commentary should begin with the # symbol " << endl;
      }
    }
  }

  //Closing the input file
  getdata.close();

  //Check that the input file contains the expected number of data.
  if (it < NXCELL * NYCELL)
  {
    cerr << rif_namefile << ": ERROR: the number of data in this file is too small!" << endl;
    exit(EXIT_FAILURE);
  }

  //Final check: Does all the grid cells were filled with a value?
  for (int i = 1; i < NXCELL + 1; i++)
  {
    for (int j = 1; j < NYCELL + 1; j++)
    {
      if (rain_c[i][j] >= MAX_SCAL)
      {
        cerr << rif_namefile << ": ERROR: the value for the point x =" << (i - 0.5) * DX << " y = " << (j - 0.5) * DY << " is missing!" << endl;
        exit(EXIT_FAILURE);
      }
    } //end for j
  }   //end for i
}

Rif_read::~Rif_read()
{
}
