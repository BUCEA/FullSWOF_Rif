/**
 * @file rain_read.cpp 
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


#include "rain_read.hpp"

Rain_read::Rain_read(Parameters & par):Rain(par){
  
  /**
   * @details Defines the name of the file for the initialization and creates two tables (times and intensity) from the data read in the file.
   * @param[in] par parameter, contains all the values from the parameters file.
   * @warning (rain_namefile): ERROR: cannot open the rain file.
   * @warning (rain_namefile): ERROR: line ***.
   * @warning (rain_namefile): ERROR: the first time must be t = 0.
   */
  
  rain_namefile = par.get_rainNameFile();
  num_lin=0;
  time_value_prev=0.;
  
  ifstream prain(rain_namefile.c_str(),ios::in);
  if (!prain){
    //The input file cannot be opened.
    cout << rain_namefile << ": ERROR: cannot open the rain file\n";
    exit(EXIT_FAILURE);
  }
  
  // As long as the end of the file is not reached, the next line is read.
  while (!prain.eof()){
    num_lin++;
    getline (prain, line, '\n');       // read a line
    istringstream entree (line);
    if (entree >> time_value >> intensity_value) {
      
      //Error if the times are decreasing.
      if(time_value < time_value_prev){
        cerr<< rain_namefile <<": ERROR: line "<< num_lin<<", t = "<< time_value <<" must be greater than " << time_value_prev <<"."<< endl;
        exit(EXIT_FAILURE);
      }
      
      //Error if the rain intensity is negative.
      if(intensity_value < 0){
        cerr<< rain_namefile <<": ERROR: line "<< num_lin<<", rain intensity = "<< intensity_value <<" must be positive."<< endl;
        exit(EXIT_FAILURE);
      }
      
      rain_times.insert(rain_times.end(),time_value);
      rain_tab.insert(rain_tab.end(),intensity_value);
      
      time_value_prev = time_value;
    }
    else{
      car='#';//Initialization of the character used to identify the beginning of a comment
      istringstream entree_car (line);
      entree_car >> car;
      if (car != '#'){
        cout << rain_namefile << ": WARNING: line "<< num_lin <<"; a commentary should begin with the # symbol " <<endl;
      }
    }
  }
  
  //Closing the input file
  prain.close();
  
  //The last time is the Maximum finite representable floating-point number,
  //so it's not necessary to write the final time in the rain file.
  rain_times.insert(rain_times.end(),DBL_MAX);
  
  if (rain_times[0] > 0.0){
    cerr<< rain_namefile <<": ERROR: the first time must be t = 0." << endl;
    exit(EXIT_FAILURE);
  }
  
  /*The index of the table rain_times is initialized to zero*/
  ind =0;
}


void Rain_read:: rain_func(SCALAR time, TAB & Tab_rain){
  
  /**
   * @details Initializes the rain to the values read in the corresponding file.
   * @param[in] time current time.
   * @param[in, out] Tab_rain rain intensity at the current time on each cell.
   * @note As the times read in the file must start with t = 0, Tab_rain is initialized.
   */
  
  if (time >= rain_times[ind]){
    for (int i=1 ; i<NXCELL+1 ; i++){
      for (int j=1 ; j<NYCELL+1 ; j++){
        Tab_rain[i][j]=(rain_tab[ind]);
      } //end for j
    } //end for i
    ind++;
  }
}


Rain_read::~Rain_read(){
}
