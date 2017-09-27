/**
 * @file parameters.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2011-2015)
 * @author Frederic Darboux <frederic.darboux@orleans.inra.fr> (2014)
 * @version 1.06.01
 * @date 2016-01-04
 *
 * @brief Gets parameters
 * @details 
 * Reads the parameters, checks their values.
 *
 * @copyright License Cecill-V2  \n
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

void Parameters::setparameters(const char *FILENAME)
{

  /**
   * @details
   * Gets all the parameters from the file FILENAME, check and affect them.
   * The values used by FullSWOF_2D are saved in the file parameters.dat. 
   * These values are also printed in the terminal when the code is run. 
   * @param[in] FILENAME name of the paramters file.
   * @warning parameters.txt: ERROR: ***.
   * @warning parameters.txt: WARNING: ***.
   * @warning ERROR: the *** file *** does not exists in the directory Inputs.
   * @warning Impossible to open the *** file. Verify if the directory *** exists.
   * @note If a value cannot be affected correctly, the code will exit with failure termination code.\n
   * If parameters.dat cannot be opened, the code will exit with failure termination code.
   */

  string path_input_directory("./Inputs/");
  string path_output_directory("./Outputs");
  Parser fileParser(FILENAME);

  /*-----------------------------Nxcell:--------------------------------------------------*/
  Nxcell = atoi(fileParser.getValue("<Nxcell>").c_str()); // atoi--convert string to integer
  if (Nxcell < 1)
  {
    cerr << " parameters.txt: ERROR: the number of cells must be greater or equal to 1." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------Nycell:--------------------------------------------------*/
  Nycell = atoi(fileParser.getValue("<Nycell>").c_str());
  if (Nycell < 1)
  {
    cerr << " parameters.txt: ERROR: the number of cells must be greater or equal to 1." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------T:-------------------------------------------------------*/
  T = atof(fileParser.getValue("<T>").c_str()); //atof--convert string to double
  if (T < 0.0)
  {
    cerr << " parameters.txt: ERROR: the final time must be non-negative." << endl;
    exit(EXIT_FAILURE);
  }
  /*-----------------------------nbtimes:-------------------------------------------------*/
  nbtimes = atoi(fileParser.getValue("<nbtimes>").c_str());
  if (nbtimes < 0)
  {
    cerr << " parameters.txt: ERROR: the number of times saved must be non-negative." << endl;
    exit(EXIT_FAILURE);
  }
  if (1 == nbtimes)
  {
    cerr << " parameters.txt: WARNING: the initial and the final time are already saved. Using nbtimes=0." << endl;
    nbtimes = 0;
  }

  /*-----------------------------scheme_type:---------------------------------------------*/
  scheme_type = atoi(fileParser.getValue("<scheme_type>").c_str());
  if (scheme_type < 1 || scheme_type > 2)
  {
    cerr << " parameters.txt: ERROR:  the choice " << scheme_type << " for the type of scheme does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------cfl_fix:-------------------------------------------------*/
  cfl_fix = atof(fileParser.getValue("<cflfix>").c_str());
  if (cfl_fix > 1.0)
  {
    cerr << " parameters.txt: WARNING: you are running the code with CFL = " << cfl_fix << endl;
  }
  if (cfl_fix < 0.0)
  {
    cerr << " parameters.txt: ERROR: the CFL must be non-negative." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------dt_fix:--------------------------------------------------*/
  dt_fix = atof(fileParser.getValue("<dtfix>").c_str());
  if (2 == scheme_type)
  {
    if (dt_fix <= 0.0)
    {
      cerr << " parameters.txt: ERROR: the time step must be non-negative." << endl;
      exit(EXIT_FAILURE);
    }
  }

  /*-----------------------------L:-------------------------------------------------------*/
  L = atof(fileParser.getValue("<L>").c_str());
  if (L <= 0.0)
  {
    cerr << " parameters.txt: ERROR: the length of the domain must be (strictly) positive." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------l:-------------------------------------------------------*/
  l = atof(fileParser.getValue("<l>").c_str());
  if (L <= 0.0)
  {
    cerr << " parameters.txt: ERROR: the width of the domain must be (strictly) positive." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------Lbound:--------------------------------------------------*/
  // type of left boundary condition
  parser_output = fileParser.getValue("<Lbound>").c_str();
  if (parser_output.size() < 1)
  {
    cerr << " parameters.txt: ERROR: the left boundary condition must be specified." << endl;
    exit(EXIT_FAILURE);
  }
  Lbound = atoi(parser_output.c_str());
  if (Lbound < 1 || Lbound > 5)
  {
    cerr << " parameters.txt: ERROR: the left boundary condition " << Lbound << " does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------Rbound:--------------------------------------------------*/
  // type of right boundary condition
  parser_output = fileParser.getValue("<Rbound>").c_str();
  if (parser_output.size() < 1)
  {
    cerr << " parameters.txt: ERROR: the right boundary condition must be specified." << endl;
    exit(EXIT_FAILURE);
  }
  Rbound = atoi(parser_output.c_str());
  if (Rbound < 1 || Rbound > 5)
  {
    cerr << " parameters.txt: ERROR: the right boundary condition " << Rbound << " does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------left_imp_discharge:--------------------------------------*/
  // left boundary condition : imposed h, q
  parser_Limp_q = fileParser.getValue("<left_imp_discharge>").c_str();
  left_imp_discharge = atof(parser_Limp_q.c_str());

  /*-----------------------------left_imp_h:----------------------------------------------*/
  parser_Limp_h = fileParser.getValue("<left_imp_h>").c_str();
  left_imp_h = atof(parser_Limp_h.c_str());

  switch (Lbound)
  { // LEFT BOUNDARY CONDITION

  case 1: // imposed water height
    if (parser_Limp_q.size() < 1)
    { // if imp_q is empty
      cerr << " parameters.txt: WARNING: imposed discharge not specified in the left boundary condition." << endl;
    }
    if (left_imp_h < 0)
    { // if negative imp_h
      cerr << " parameters.txt: ERROR: negative imposed water height in the left boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;

  // case 2: // wall condition
  // nothing to test on these values

  // case 3: // Neumann condition
  // nothing to test on these values

  case 4: // periodic case
    if (Rbound != 4)
    { // if the other bound is not periodic
      cerr << " parameters.txt: ERROR: you must choose a periodic right condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;

  case 5: // imposed discharge
    if (parser_Limp_q.size() < 1)
    { // if imp_q is empty
      cerr << " parameters.txt: ERROR: imposed discharge not specified in the left boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    if (parser_Limp_h.size() < 1)
    { // if imp_h is empty
      cerr << " parameters.txt: WARNING: imposed water height not specified in the left boundary condition." << endl;
    }
    if (left_imp_h < 0)
    { // if negative imp_h
      cerr << " parameters.txt: ERROR: negative imposed water height in the left boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;
  }

  /*-----------------------------right_imp_discharge:-------------------------------------*/
  parser_Rimp_q = fileParser.getValue("<right_imp_discharge>").c_str();
  right_imp_discharge = atof(parser_Rimp_q.c_str());

  /*-----------------------------right_imp_h:---------------------------------------------*/
  parser_Rimp_h = fileParser.getValue("<right_imp_h>").c_str();
  right_imp_h = atof(parser_Rimp_h.c_str());

  switch (Rbound)
  { // RIGHT BOUNDARY CONDITION

  case 1: // imposed water height
    if (parser_Rimp_q.size() < 1)
    { // if imp_q is empty
      cerr << " parameters.txt: WARNING: imposed discharge not specified in the right boundary condition." << endl;
    }
    if (right_imp_h < 0)
    { // if negative imp_h
      cerr << " parameters.txt: ERROR: negative imposed water height in the right boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;

  // case 2: // wall condition
  // nothing to test on these values

  // case 3: // Neumann condition
  // nothing to test on these values

  case 4: // periodic case
    if (Lbound != 4)
    { // if the other bound is not periodic
      cerr << " parameters.txt: ERROR: you must choose a periodic left condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;

  case 5: // imposed discharge
    if (parser_Rimp_q.size() < 1)
    { // if imp_q is empty
      cerr << " parameters.txt: ERROR: imposed discharge not specified in the right boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    if (parser_Rimp_h.size() < 1)
    { // if imp_h is empty
      cerr << " parameters.txt: WARNING: imposed water height not specified in the right boundary condition." << endl;
    }
    if (right_imp_h < 0)
    { // if negative imp_h
      cerr << " parameters.txt: ERROR: negative imposed water height in the right boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;
  }

  /*-----------------------------Bbound:--------------------------------------------------*/
  // type of bottom boundary condition
  parser_output = fileParser.getValue("<Bbound>").c_str();
  if (parser_output.size() < 1)
  {
    cerr << " parameters.txt: ERROR: the bottom boundary condition must be specified." << endl;
    exit(EXIT_FAILURE);
  }
  Bbound = atoi(parser_output.c_str());
  if (Bbound < 1 || Bbound > 5)
  {
    cerr << " parameters.txt: ERROR: the bottom boundary condition " << Bbound << " does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------Tbound:--------------------------------------------------*/
  // type of top boundary condition
  parser_output = fileParser.getValue("<Tbound>").c_str();
  if (parser_output.size() < 1)
  {
    cerr << " parameters.txt: ERROR: the top boundary condition must be specified." << endl;
    exit(EXIT_FAILURE);
  }
  Tbound = atoi(parser_output.c_str());
  if (Tbound < 1 || Tbound > 5)
  {
    cerr << " parameters.txt: ERROR: the bottom boundary condition " << Tbound << " does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------bottom_imp_discharge:------------------------------------*/
  parser_Bimp_q = fileParser.getValue("<bottom_imp_discharge>").c_str();
  bottom_imp_discharge = atof(parser_Bimp_q.c_str());

  /*-----------------------------bottom_imp_h:--------------------------------------------*/
  parser_Bimp_h = fileParser.getValue("<bottom_imp_h>").c_str();
  bottom_imp_h = atof(parser_Bimp_h.c_str());

  switch (Bbound)
  { // BOTTOM BOUNDARY CONDITION

  case 1: // imposed water height
    if (parser_Bimp_q.size() < 1)
    { // if imp_q is empty
      cerr << " parameters.txt: WARNING: imposed discharge not specified in the bottom boundary condition." << endl;
    }
    if (bottom_imp_h < 0)
    { // if negative imp_h
      cerr << " parameters.txt: ERROR: negative imposed water height in the bottom boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;

  // case 2: // wall condition
  // nothing to test on these values

  // case 3: // Neumann condition
  // nothing to test on these values

  case 4: // periodic case
    if (Tbound != 4)
    { // if the other bound is not periodic
      cerr << " parameters.txt: ERROR: you must choose a periodic top condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;

  case 5: // imposed discharge
    if (parser_Bimp_q.size() < 1)
    { // if imp_q is empty
      cerr << " parameters.txt: ERROR: imposed discharge not specified in the bottom boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    if (parser_Bimp_h.size() < 1)
    { // if imp_h is empty
      cerr << " parameters.txt: WARNING: imposed water height not specified in the bottom boundary condition." << endl;
    }
    if (bottom_imp_h < 0)
    { // if negative imp_h
      cerr << " parameters.txt: ERROR: negative imposed water height in the bottom boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;
  }

  /*-----------------------------top_imp_discharge:---------------------------------------*/
  parser_Timp_q = fileParser.getValue("<top_imp_discharge>").c_str();
  top_imp_discharge = atof(parser_Timp_q.c_str());

  /*-----------------------------top_imp_h:-----------------------------------------------*/
  parser_Timp_h = fileParser.getValue("<top_imp_h>").c_str();
  top_imp_h = atof(parser_Timp_h.c_str());

  switch (Tbound)
  { // TOP BOUNDARY CONDITION

  case 1: // imposed water height
    if (parser_Timp_q.size() < 1)
    { // if imp_q is empty
      cerr << " parameters.txt: WARNING: imposed discharge not specified in the top boundary condition." << endl;
    }
    if (top_imp_h < 0)
    { // if negative imp_h
      cerr << " parameters.txt: ERROR: negative imposed water height in the top boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;

  // case 2: // wall condition
  // nothing to test on these values

  // case 3: // Neumann condition
  // nothing to test on these values

  case 4: // periodic case
    if (Bbound != 4)
    { // if the other bound is not periodic
      cerr << " parameters.txt: ERROR: you must choose a periodic bottom condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;

  case 5: // imposed discharge
    if (parser_Timp_q.size() < 1)
    { // if imp_q is empty
      cerr << " parameters.txt: ERROR: imposed discharge not specified in the top boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    if (parser_Timp_h.size() < 1)
    { // if imp_h is empty
      cerr << " parameters.txt: WARNING: imposed water height not specified in the top boundary condition." << endl;
    }
    if (top_imp_h < 0)
    { // if negative imp_h
      cerr << " parameters.txt: ERROR: negative imposed water height in the top boundary condition." << endl;
      exit(EXIT_FAILURE);
    }
    break;
  }

  /*-----------------------------fric:----------------------------------------------------*/
  parser_output = fileParser.getValue("<fric>").c_str();
  if (parser_output.size() < 1)
  { // if fric is empty
    cerr << " parameters.txt: ERROR: the friction law  must be specified." << endl;
    exit(EXIT_FAILURE);
  }

  fric = atoi(parser_output.c_str());

  if (fric < 0 || fric > 3)
  {
    cerr << " parameters.txt: ERROR: the friction law " << fric << " does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  if (fric > 0)
  { //if the case differs from no friction
    if (fileParser.getValue("<fric_init>").size() < 1)
    { // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of friction term (<fric_init>)." << endl;
      exit(EXIT_FAILURE);
    }

    fric_init = atoi(fileParser.getValue("<fric_init>").c_str());

    if (fric_init < 1 || fric_init > 2)
    {
      cerr << " parameters.txt: ERROR: initialization number " << fric_init << " for friction term (<fric_init>) does not exist." << endl;
      exit(EXIT_FAILURE);
    }

    if (1 == fric_init)
    { //initialization from a file
      fric_NF = fileParser.getValue("<fric_NF>");
      if (fric_NF.size() < 1)
      { // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of friction." << endl;
        exit(EXIT_FAILURE);
      }
      fric_namefile = path_input_directory + fric_NF;
      if (-1 == access(fric_namefile.c_str(), R_OK))
      {
        cerr << " parameters.txt: ERROR: the " << fric_NF << " file does not exists in the directory Inputs." << endl;
        exit(EXIT_FAILURE);
      }
    }
    else
    {

      /*-----------------------------fric_coef:------------------------------------------------*/
      friccoef = atof(fileParser.getValue("<friccoef>").c_str());

      if (friccoef < 0)
      { //if negative
        cerr << " parameters.txt: ERROR: negative friction coefficient." << endl;
        exit(EXIT_FAILURE);
      }
      if (0. >= friccoef)
      {
        cerr << " parameters.txt: WARNING: friction law with null friction coefficient " << endl;
        cerr << "*********************************************************" << endl;
      }
    }
  }

  /*-----------------------------flux:----------------------------------------------------*/
  parser_output = fileParser.getValue("<flux>").c_str();
  if (parser_output.size() < 1)
  { // if flux is empty
    cerr << " parameters.txt: WARNING: the numerical flux is empty, using HLL instead" << endl;
    cerr << "*********************************************************" << endl;
    flux = 2;
  }
  else
  {
    flux = atoi(parser_output.c_str());
    if (flux < 1 || flux > 5)
    {
      cerr << " parameters.txt: ERROR: the numerical flux " << flux << " does not exist." << endl;
      exit(EXIT_FAILURE);
    }
  }

  /*-----------------------------order:---------------------------------------------------*/
  parser_output = fileParser.getValue("<order>").c_str();
  if (parser_output.size() < 1)
  { // if order is empty
    cerr << " parameters.txt: WARNING: the order of the scheme is empty, using order 2 instead" << endl;
    cerr << "*********************************************************" << endl;
    order = 2;
  }
  else
  {
    order = atoi(parser_output.c_str());

    if (order < 1 || order > 2)
    {
      cerr << "parameters.txt: ERROR: the order " << order << " does not exist for this scheme." << endl;
      exit(EXIT_FAILURE);
    }
  }

  /*-----------------------------rec:-----------------------------------------------------*/
  if (2 == order)
  { // in case of order 2 only
    parser_output = fileParser.getValue("<rec>").c_str();
    if (parser_output.size() < 1)
    { // if rec is empty
      cerr << " parameters.txt: WARNING: the reconstruction is empty, using MUSCL instead" << endl;
      cerr << "*********************************************************" << endl;
      rec = 1; //using MUSCL by default
    }
    else
    {
      rec = atoi(parser_output.c_str());
      if (rec < 1 || rec > 3)
      {
        cerr << " parameters.txt: ERROR: the reconstruction " << rec << " does not exist." << endl;
        exit(EXIT_FAILURE);
      }

      if (rec > 1)
      { //in case of ENO or ENOmod
        /*-----------------------------amortENO:------------------------------------------------*/
        parser_output2 = fileParser.getValue("<amortENO>").c_str();
        if (parser_output2.size() < 1)
        { // if amortENO is empty
          cerr << " parameters.txt: ERROR: the AmortENO is empty" << endl;
          cerr << "*********************************************************" << endl;
          exit(EXIT_FAILURE);
        }
        else
        {
          amortENO = atof(parser_output2.c_str());
          if (amortENO < 0 || amortENO > 1)
          {
            cerr << "parameters.txt: ERROR: the AmortENO " << amortENO << " is not between 0 and 1." << endl;
            exit(EXIT_FAILURE);
          }
        }
        if (3 == rec)
        { //in case of  ENOmod
          /*-----------------------------modifENO:------------------------------------------------*/
          parser_output2 = fileParser.getValue("<modifENO>").c_str();
          if (parser_output2.size() < 1)
          { // if modifENO is empty
            cerr << " parameters.txt: ERROR: the modifENO is empty" << endl;
            cerr << "*********************************************************" << endl;
            exit(EXIT_FAILURE);
          }
          else
          {
            modifENO = atof(parser_output2.c_str());
            if (modifENO < 0 || modifENO > 1)
            {
              cerr << " parameters.txt: ERROR: the modifENO " << modifENO << " is not between 0 and 1." << endl;
              exit(EXIT_FAILURE);
            }
          }
        } //endif rec==3

      } //endif in case of ENO or ENOmod
    }   //end rec treatment

    /*-----------------------------lim:-----------------------------------------------------*/
    parser_output = fileParser.getValue("<lim>").c_str();
    if (parser_output.size() < 1)
    { // if lim is empty
      cerr << " parameters.txt: WARNING: the limiter is empty, using Minmod instead" << endl;
      cerr << "*********************************************************" << endl;
      lim = 1; //using Minmod by default
    }
    else
    {
      lim = atoi(parser_output.c_str());
      if (lim < 1 || lim > 3)
      {
        cerr << "parameters.txt: ERROR: the limiter " << lim << " does not exist." << endl;
        exit(EXIT_FAILURE);
      }
    }
  } //endif order 2

  /*-----------------------------inf:-----------------------------------------------------*/
  parser_output = fileParser.getValue("<inf>").c_str();
  if (parser_output.size() < 1)
  { // if inf is empty
    cerr << " parameters.txt: ERROR: the infiltration model  must be specified." << endl;
    exit(EXIT_FAILURE);
  }

  inf = atoi(parser_output.c_str());

  if (inf < 0 || inf > 1)
  {
    cerr << " parameters.txt: ERROR: the infiltration model " << inf << " does not exist." << endl;
    exit(EXIT_FAILURE);
  }
  if (inf > 0)
  { //if the case differs from no infiltration
    /*-----------------------------Kc:------------------------------------------------------*/
    if (fileParser.getValue("<Kc_init>").size() < 1)
    { // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of hydraulic conductivity of the crust (<Kc_init>)." << endl;
      exit(EXIT_FAILURE);
    }

    Kc_init = atoi(fileParser.getValue("<Kc_init>").c_str());

    if (Kc_init < 1 || Kc_init > 2)
    {
      cerr << " parameters.txt: ERROR: initialization number " << Kc_init << " for hydraulic conductivity of the crust (<Kc_init>) does not exist." << endl;
      exit(EXIT_FAILURE);
    }

    if (1 == Kc_init)
    { //initialization from a file
      Kc_NF = fileParser.getValue("<Kc_NF>");
      if (Kc_NF.size() < 1)
      { // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of Kc." << endl;
        exit(EXIT_FAILURE);
      }
      Kc_namefile = path_input_directory + Kc_NF;
      if (-1 == access(Kc_namefile.c_str(), R_OK))
      {
        cerr << " parameters.txt: ERROR: the " << Kc_NF << " file does not exists in the directory Inputs." << endl;
        exit(EXIT_FAILURE);
      }
    }
    else
    {
      Kc_coef = atof(fileParser.getValue("<Kccoef>").c_str());

      if (Kc_coef < 0)
      { //if hydraulic conductivity of the crust negative
        cerr << " parameters.txt: ERROR: negative hydraulic conductivity of the crust." << endl;
        exit(EXIT_FAILURE);
      }

      if (0. >= Kc_coef)
      {
        cerr << " parameters.txt: WARNING: infiltration model with null hydraulic conductivity of the crust  " << endl;
        cerr << "*********************************************************" << endl;
      }
    }
    /*-----------------------------Ks:------------------------------------------------------*/
    if (fileParser.getValue("<Ks_init>").size() < 1)
    { // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of hydraulic conductivity of the soil (<Ks_init>)." << endl;
      exit(EXIT_FAILURE);
    }

    Ks_init = atoi(fileParser.getValue("<Ks_init>").c_str());

    if (Ks_init < 1 || Ks_init > 2)
    {
      cerr << " parameters.txt: ERROR: initialization number " << Ks_init << " for hydraulic conductivity of the soil (<Ks_init>) does not exist." << endl;
      exit(EXIT_FAILURE);
    }
    if (1 == Ks_init)
    { //initialization from a file
      Ks_NF = fileParser.getValue("<Ks_NF>");
      if (Ks_NF.size() < 1)
      { // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of Ks." << endl;
        exit(EXIT_FAILURE);
      }
      Ks_namefile = path_input_directory + Ks_NF;
      if (-1 == access(Ks_namefile.c_str(), R_OK))
      {
        cerr << " parameters.txt: ERROR: the " << Ks_NF << " file does not exists in the directory Inputs." << endl;
        exit(EXIT_FAILURE);
      }
    }
    else
    {
      Ks_coef = atof(fileParser.getValue("<Kscoef>").c_str());
      if (Ks_coef < 0)
      { //if hydraulic conductivity of the soil negative
        cerr << " parameters.txt: ERROR: negative hydraulic conductivity of the soil." << endl;
        exit(EXIT_FAILURE);
      }
    }

    /*-----------------------------dtheta:--------------------------------------------------*/
    if (fileParser.getValue("<dtheta_init>").size() < 1)
    { // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of water content (<dtheta_init>)." << endl;
      exit(EXIT_FAILURE);
    }

    dtheta_init = atoi(fileParser.getValue("<dtheta_init>").c_str());

    if (dtheta_init < 1 || dtheta_init > 2)
    {
      cerr << " parameters.txt: ERROR: initialization number " << dtheta_init << " for water content (<dtheta_init>) does not exist." << endl;
      exit(EXIT_FAILURE);
    }
    if (1 == dtheta_init)
    { //initialization from a file
      dtheta_NF = fileParser.getValue("<dtheta_NF>");
      if (dtheta_NF.size() < 1)
      { // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of dtheta." << endl;
        exit(EXIT_FAILURE);
      }
      dtheta_namefile = path_input_directory + dtheta_NF;
      if (-1 == access(dtheta_namefile.c_str(), R_OK))
      {
        cerr << " parameters.txt: ERROR: the " << dtheta_NF << " file does not exists in the directory Inputs." << endl;
        exit(EXIT_FAILURE);
      }
    }
    else
    {
      dtheta_coef = atof(fileParser.getValue("<dthetacoef>").c_str());
      if (dtheta_coef < 0)
      { //if water content negative
        cerr << " parameters.txt: ERROR: negative water content." << endl;
        exit(EXIT_FAILURE);
      }

      if (dtheta_coef > 1)
      { //if water content are above 1
        cerr << " parameters.txt: WARNING: the value of dtheta seems very large." << endl;
      }
    }

    /*-----------------------------Psi:-----------------------------------------------------*/
    if (fileParser.getValue("<Psi_init>").size() < 1)
    { // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of load pressure (<Psi_init>)." << endl;
      exit(EXIT_FAILURE);
    }

    Psi_init = atoi(fileParser.getValue("<Psi_init>").c_str());

    if (Psi_init < 1 || Psi_init > 2)
    {
      cerr << " parameters.txt: ERROR: initialization number " << Psi_init << " for load pressure (<Psi_init>) does not exist." << endl;
      exit(EXIT_FAILURE);
    }
    if (1 == Psi_init)
    { //initialization from a file
      Psi_NF = fileParser.getValue("<Psi_NF>");
      if (Psi_NF.size() < 1)
      { // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of Psi." << endl;
        exit(EXIT_FAILURE);
      }
      Psi_namefile = path_input_directory + Psi_NF;
      if (-1 == access(Psi_namefile.c_str(), R_OK))
      {
        cerr << " parameters.txt: ERROR: the " << Psi_NF << " file does not exists in the directory Inputs." << endl;
        exit(EXIT_FAILURE);
      }
    }
    else
    {
      Psi_coef = atof(fileParser.getValue("<Psicoef>").c_str());
      if (Psi_coef < 0)
      { //if load pressure negative
        cerr << " parameters.txt: ERROR: negative load pressure." << endl;
        exit(EXIT_FAILURE);
      }
    }

    /*-----------------------------zcrust:--------------------------------------------------*/
    if (fileParser.getValue("<zcrust_init>").size() < 1)
    { // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of thickness of the crust (<zcrust_init>)." << endl;
      exit(EXIT_FAILURE);
    }

    zcrust_init = atoi(fileParser.getValue("<zcrust_init>").c_str());

    if (zcrust_init < 1 || zcrust_init > 2)
    {
      cerr << " parameters.txt: ERROR: initialization number " << zcrust_init << " for thickness of the crust (<zcrust_init>) does not exist." << endl;
      exit(EXIT_FAILURE);
    }
    if (1 == zcrust_init)
    { //initialization from a file
      zcrust_NF = fileParser.getValue("<zcrust_NF>");
      if (zcrust_NF.size() < 1)
      { // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of zcrust." << endl;
        exit(EXIT_FAILURE);
      }
      zcrust_namefile = path_input_directory + zcrust_NF;
      if (-1 == access(zcrust_namefile.c_str(), R_OK))
      {
        cerr << " parameters.txt: ERROR: the " << zcrust_NF << " file does not exists in the directory Inputs." << endl;
        exit(EXIT_FAILURE);
      }
    }
    else
    {
      zcrust_coef = atof(fileParser.getValue("<zcrustcoef>").c_str());
      if (zcrust_coef < 0)
      { //if thickness of the crust negative
        cerr << " parameters.txt: ERROR: negative thickness of the crust." << endl;
        exit(EXIT_FAILURE);
      }
    }

    /*-----------------------------imax:----------------------------------------------------*/
    if (fileParser.getValue("<imax_init>").size() < 1)
    { // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of Maximun infiltration rate (<imax_init>)." << endl;
      exit(EXIT_FAILURE);
    }

    imax_init = atoi(fileParser.getValue("<imax_init>").c_str());

    if (imax_init < 1 || imax_init > 2)
    {
      cerr << " parameters.txt: ERROR: initialization number " << imax_init << " for Maximun infiltration rate (<imax_init>) does not exist." << endl;
      exit(EXIT_FAILURE);
    }
    if (1 == imax_init)
    { //initialization from a file
      imax_NF = fileParser.getValue("<imax_NF>");
      if (imax_NF.size() < 1)
      { // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to have a file for this choice of initialization of imax." << endl;
        exit(EXIT_FAILURE);
      }
      imax_namefile = path_input_directory + imax_NF;
      if (-1 == access(imax_namefile.c_str(), R_OK))
      {
        cerr << " parameters.txt: ERROR: the " << imax_NF << " file does not exists in the directory Inputs." << endl;
        exit(EXIT_FAILURE);
      }
    }
    else
    {
      imax_coef = atof(fileParser.getValue("<imaxcoef>").c_str());
      if (imax_coef < 0)
      { //if Maximun infiltration rate negative
        cerr << " parameters.txt: ERROR: negative Maximun infiltration rate." << endl;
        exit(EXIT_FAILURE);
      }
    }
  } //endif the case differs from no infiltration

  /*-----------------------------topo:----------------------------------------------------*/
  topo = atoi(fileParser.getValue("<topo>").c_str());
  if (topo < 1 || topo > 3)
  {
    cerr << " parameters.txt: ERROR: initialization number " << topo << " for the topography does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------topo_NF:-------------------------------------------------*/
  topo_NF = fileParser.getValue("<topo_NF>");
  topography_namefile = path_input_directory + topo_NF;
  if (1 == topo && -1 == access(topography_namefile.c_str(), R_OK))
  {
    cerr << " parameters.txt: ERROR: the topography file " << topo_NF << " does not exists in the directory Inputs." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------huv_init:------------------------------------------------*/
  huv_init = atoi(fileParser.getValue("<huv_init>").c_str());
  if (huv_init < 1 || huv_init > 5)
  {
    cerr << " parameters.txt: ERROR: initialization number " << huv_init << " for h, u does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------huv_NF:--------------------------------------------------*/
  huv_NF = fileParser.getValue("<huv_NF>");
  huv_namefile = path_input_directory + huv_NF;
  if (1 == huv_init && -1 == access(huv_namefile.c_str(), R_OK))
  {
    cerr << " parameters.txt: ERROR: the huv file " << huv_NF << " does not exists in the directory Inputs." << endl;
    exit(EXIT_FAILURE);
  }

  /*------------------------------rif_init:-----------------------------------------------*/
  rif_init = atoi(fileParser.getValue("<rif_init>").c_str());
  if (rif_init <1 || rif_init >5)
  {
    cerr << " parameters.txt: Error: initialization number " << rif_init << " for r, i and f does not exist." <<endl;
    exit(EXIT_FAILURE);
  }

  /*------------------------------rif_NF:-------------------------------------------------*/
  rif_NF = fileParser.getValue("<rif_NF>");
  rif_namefile = path_input_directory + rif_NF;
  if (1 == rif_init && -1 == access(rif_namefile.c_str(), R_OK))
  {
    cerr << " parameters.txt: Error: the rif file " << rif_NF << " does not exists in the directory Inputs."<< endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------rain:----------------------------------------------------*/
  rain = atoi(fileParser.getValue("<rain>").c_str());
  if (rain < 0 || rain > 2)
  {
    cerr << " parameters.txt: ERROR: the choice " << rain << " for rain does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------rain_NF:-------------------------------------------------*/
  rain_NF = fileParser.getValue("<rain_NF>");
  rain_namefile = path_input_directory + rain_NF;
  if (1 == rain && -1 == access(rain_namefile.c_str(), R_OK))
  {
    cerr << " parameters.txt: ERROR: the rain file " << rain_NF << " does not exists in the directory Inputs" << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------suffix_o:------------------------------------------------*/
  suffix_outputs = fileParser.getValue("<suffix_o>");

  /*-----------------------------output_f:--------------------------------------------------------*/
  if (0 == nbtimes)
  {
    output_format = 0;
  }
  else
  {
    output_format = atoi(fileParser.getValue("<output_f>").c_str());
    if (output_format < 1 || output_format > 2)
    {
      cerr << " parameters.txt: ERROR: the choice " << output_format << " for the format of the output file does not exist." << endl;
      exit(EXIT_FAILURE);
    }
  }

  dx = L / Nxcell; //dx and dy calculated
  dy = l / Nycell;

  output_directory = path_output_directory + suffix_outputs + "/";

  namefile = output_directory + "parameters.dat";
  ofstream param(namefile.c_str(), ios::out);        //ios::out open for output operation
  if (!param)
  {
    cerr << " parameters.txt: ERROR: Impossible to open the " << namefile.c_str() << " file\n";
    cerr << " parameters.txt: ERROR: Verify if the directory " << output_directory << " exists\n";
    exit(EXIT_FAILURE);
  }

  param << "#####################################################################" << endl;
  param << "# Parameters used in FullSWOF_2D software." << endl;
  param << "#####################################################################" << endl;
  param << endl;
  param << "Number of meshes (x-axis)  <Nxcell>:: " << Nxcell << endl;
  param << "Number of meshes (y-axis)  <Nycell>:: " << Nycell << endl;
  param << endl;
  param << "Time of simulation  <T>:: " << T << endl;
  param << "Number of times saved <nbtimes>:: " << nbtimes << endl;
  param << endl;
  param << "Choice of type of scheme (1=fixed cfl  2=fixed dt) <scheme_type>:: " << scheme_type << endl;
  if (2 == scheme_type)
  { // if fixed dt
    param << "Timestep (in seconds) <dtfix>:: " << dt_fix << endl;
  }
  else
  {
    param << "Timestep (in seconds) <dtfix>:: " << endl;
  }
  param << "Value of the cfl  <cflfix>:: " << cfl_fix << endl;
  param << endl;
  param << "Length of the domain in respect to x  <L>:: " << L << endl;
  param << "Length of the domain in respect to y  <l>:: " << l << endl;
  param << endl;
  param << endl;
  param << "Left Boundary condition   (x = 0)   (1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q)  <Lbound>:: " << Lbound << endl;
  if (2 == Lbound || 3 == Lbound || 4 == Lbound)
  {                                                                           // if wall, Neumann,  or periodic
    param << "Imposed discharge in left bc <left_imp_discharge> :: " << endl; // empty as nothing is used
    param << "Imposed water height in left bc <left_imp_h> :: " << endl;      // empty as nothing is used
  }
  else
  {
    param << "Imposed discharge in left bc <left_imp_discharge> :: " << left_imp_discharge << endl;     // writes 0 if the value was empty
    param << "Imposed height in left bc (if flow supercritical) <left_imp_h>:: " << left_imp_h << endl; // writes 0 if the value was empty
  }
  param << endl;
  param << "Right Boundary condition  (x = xmax)(1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q)  <Rbound>:: " << Rbound << endl;
  if (2 == Rbound || 3 == Rbound || 4 == Rbound)
  {                                                                             // if wall, Neumann,  or periodic
    param << "Imposed discharge in right bc <right_imp_discharge> :: " << endl; // empty as nothing is used
    param << "Imposed water height in right bc <right_imp_h> :: " << endl;      // empty as nothing is used
  }
  else
  {
    param << "Imposed discharge in right bc <right_imp_discharge> :: " << right_imp_discharge << endl;     // writes 0 if the value was empty
    param << "Imposed height in right bc (if flow supercritical) <right_imp_h>:: " << right_imp_h << endl; // writes 0 if the value was empty
  }
  param << endl;
  param << "Bottom Boundary condition (y = 0)   (1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q)  <Bbound>:: " << Bbound << endl;
  if (2 == Bbound || 3 == Bbound || 4 == Bbound)
  {                                                                               // if wall, Neumann,  or periodic
    param << "Imposed discharge in bottom bc <bottom_imp_discharge> :: " << endl; // empty as nothing is used
    param << "Imposed water height in bottom bc <bottom_imp_h> :: " << endl;      // empty as nothing is used
  }
  else
  {
    param << "Imposed discharge in bottom  bc <bottom_imp_discharge> :: " << bottom_imp_discharge << endl;    // empty as nothing is used
    param << "Imposed height in bottom bc (if flow supercritical) <bottom_imp_h>:: " << bottom_imp_h << endl; // empty as nothing is used
  }
  param << endl;
  param << "Top Boundary condition (y = ymax)(1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q)  <Tbound>:: " << Tbound << endl;
  if (2 == Tbound || 3 == Tbound || 4 == Tbound)
  {                                                                         // if wall, Neumann,  or periodic
    param << "Imposed discharge in top bc <top_imp_discharge> :: " << endl; // empty as nothing is used
    param << "Imposed water height in top bc <top_imp_h> :: " << endl;      // empty as nothing is used
  }
  else
  {
    param << "Imposed discharge in top bc <top_imp_discharge> :: " << top_imp_discharge << endl;     // empty as nothing is used
    param << "Imposed height in top bc (if flow supercritical) <top_imp_h>:: " << top_imp_h << endl; // empty as nothing is used
  }
  param << endl;
  param << "Initialization of Friction (1=file 2=const_coef) <fric_init>:: " << fric_init << endl;
  param << "Friction law (0=No Friction 1=Manning 2=Darcy-Weisbach 3=laminar)  <fric>:: " << fric << endl;
  if (fric_init == 1)
  {
    param << "Name of the friction file <fric_NF>:: " << fric_NF << endl;
  }
  else
  {
    param << "Name of the friction file <fric_NF>:: " << endl;
  }
  if (0 == fric)
  {
    param << "Friction coefficient  <friccoef>:: " << endl;
  }
  else
  {
    param << "Friction coefficient  <friccoef>:: " << friccoef << endl;
  }
  param << endl;

  param << "Numerical flux (1=Rusanov 2=HLL 3=HLL2  4=HLLC 5=HLLC2)  <flux>:: " << flux << endl;
  param << endl;
  param << "Order of the scheme (1=order1 2=order2)  <order>:: " << order << endl;
  param << endl;
  if (2 == order)
  { // in case of order 2
    param << "Reconstruction (1=MUSCL 2=ENO 3=ENOmod)  <rec>:: " << rec << endl;
    if (1 == rec)
    {
      param << "AmortENO (Between 0 and 1.)  <amortENO>:: " << endl;
      param << "ModifENO (Between 0 and 1.)  <modifENO>:: " << endl;
    }
    else
    {
      param << "AmortENO (Between 0 and 1.)  <amortENO>:: " << amortENO << endl;
      if (3 == rec)
      {
        param << "ModifENO (Between 0 and 1.)  <modifENO>:: " << modifENO << endl;
      }
      else
      {
        param << "ModifENO (Between 0 and 1.)  <modifENO>:: " << endl;
      }
    }
    param << "Limiter (1=Minmod 2=VanAlbada 3=VanLeer)  <lim>:: " << lim << endl;
  }
  else
  {
    param << "Reconstruction (1=MUSCL 2=ENO 3=ENOmod)  <rec>:: " << endl;
    param << "AmortENO (Between 0 and 1.)  <amortENO>:: " << endl;
    param << "ModifENO (Between 0 and 1.)  <modifENO>:: " << endl;
    param << "Limiter (1=Minmod 2=VanAlbada 3=VanLeer)  <lim>:: " << endl;
  }

  param << endl;

  param << "Infiltration model (0=No Infiltration 1=Green-Ampt)  <inf>:: " << inf << endl;
  if (0 == inf)
  {
    param << "zcrust, thickness of the crust  (1=file 2=const_coef) <zcrust_init>:: " << endl;
    param << "zcrust coefficient <zcrustcoef>:: " << endl;
    param << "Name of the zcrust file <zcrust_NF>::" << endl;
    param << endl;

    param << "Kc, hydraulic conductivity (saturation) of the crust (1=file 2=const_coef) <Kc_init>:: " << endl;
    param << "Kc coefficient  <Kccoef>::" << endl;
    param << "Name of the Kc file <Kc_NF>::" << endl;
    param << endl;

    param << "Ks, hydraulic conductivity (saturation) of the soil (1=file 2=const_coef) <Ks_init_init>:: " << endl;
    param << "Ks coefficient  <Kscoef>::" << endl;
    param << "Name of the Ks file <Ks_NF>:: " << endl;
    param << endl;

    param << "dtheta, water content  (1=file 2=const_coef) <dtheta_init>:: " << endl;
    param << "dtheta coefficient  <dthetacoef>::" << endl;
    param << "Name of the dtheta file <dtheta_NF>::" << endl;
    param << endl;

    param << "Psi, load pressure  (1=file 2=const_coef) (1=file 2=const_coef) <Psi_init>:: " << endl;
    param << "Psi coefficient   <Psicoef>::" << endl;
    param << "Name of the dtheta file <dtheta_NF>::" << endl;
    param << endl;

    param << "imax, Maximun infiltration  rate  (1=file 2=const_coef) <imax_init>:: " << endl;
    param << "imax coefficient   <imaxcoef>::" << endl;
    param << "Name of the imax file <imax_NF>::" << endl;
  }
  else
  {
    param << endl;

    param << "zcrust, thickness of the crust (1=file 2=const_coef)  <zcrust_init>:: " << zcrust_init << endl;
    if (zcrust_init == 1)
    {
      param << "zcrust coefficient <zcrustcoef>:: " << endl;
      param << "Name of the zcrust file <zcrust_NF>:: " << zcrust_NF << endl;
    }
    else
    {
      param << "zcrust coefficient <zcrustcoef>:: " << zcrust_coef << endl;
      param << "Name of the zcrust file <zcrust_NF>:: " << endl;
    }
    param << endl;

    param << "Kc, hydraulic conductivity (saturation) of the crust (1=file 2=const_coef) <Kc_init>:: " << Kc_init << endl;
    if (Kc_init == 1)
    {
      param << "Kc coefficient  <Kccoef>:: " << endl;
      param << "Name of the Kc file <Kc_NF>::  " << Kc_NF << endl;
    }
    else
    {
      param << "Kc coefficient  <Kccoef>:: " << Ks_coef << endl;
      param << "Name of the Kc file <Kc_NF>::  " << endl;
    }
    param << endl;

    param << "Ks, hydraulic conductivity (saturation) of the soil (1=file 2=const_coef) <Ks_init>:: " << Ks_init << endl;
    if (Ks_init == 1)
    {
      param << "Ks coefficient  <Kscoef>:: " << endl;
      param << "Name of the Ks file <Ks_NF>::  " << Ks_NF << endl;
    }
    else
    {
      param << "Ks coefficient  <Kscoef>:: " << Ks_coef << endl;
      param << "Name of the Ks file <Ks_NF>::  " << endl;
    }
    param << endl;

    param << "dtheta, water content  (1=file 2=const_coef) <dtheta_init>:: " << dtheta_init << endl;
    if (dtheta_init == 1)
    {
      param << "dtheta coefficient  <dthetacoef>:: " << endl;
      param << "Name of the dtheta file <dtheta_NF>::  " << dtheta_NF << endl;
    }
    else
    {
      param << "dtheta coefficient  <dthetacoef>:: " << dtheta_coef << endl;
      param << "Name of the dtheta file <dtheta_NF>::  " << endl;
    }
    param << endl;

    param << "Psi, load pressure (1=file 2=const_coef) <Psi_init>:: " << Psi_init << endl;
    if (Psi_init == 1)
    {
      param << "Psi coefficient   <Psicoef>:: " << endl;
      param << "Name of the Psi file <Psi_NF>::  " << Psi_NF << endl;
    }
    else
    {
      param << "Psi coefficient   <Psicoef>:: " << Psi_coef << endl;
      param << "Name of the Psi file <Psi_NF>::  " << endl;
    }
    param << endl;

    param << "imax, Maximum infiltration  rate  (1=file 2=const_coef) <imax_init>:: " << imax_init << endl;
    if (imax_init == 1)
    {
      param << "imax coefficient   <imaxcoef>:: " << endl;
      param << "Name of the imax file <imax_NF>::  " << imax_NF << endl;
    }
    else
    {
      param << "imax coefficient   <imaxcoef>:: " << imax_coef << endl;
      param << "Name of the imax file <imax_NF>::  " << endl;
    }
  }
  param << endl;
  param << "Topography (1=file 2=flat 3=Thacker)  <topo>:: " << topo << endl;
  param << "Name of the topography file  <topo_NF>:: " << topo_NF << endl;
  param << endl;
  param << "Initialization of h, u and v (1=file 2=h,u&v=0 3=Thacker 4=Radial_Dam_dry 5=Radial_Dam_wet)  <huv_init>:: " << huv_init << endl;
  param << "Name of the huv initialization file  <huv_NF>:: " << huv_NF << endl;
  param << endl;
  param << "Initialization of r, i and f (1=file 2=r,i&f=0 3=Thacker 4=Radial_Dam_dry 5=Radial_Dam_wet)  <rif_init>:: " << rif_init << endl;
  param << "Name of the rif initialization file  <rif_NF>:: " << rif_NF << endl;
  param << endl;
  param << "Rain (0=no rain 1=file 2=function)  <rain>:: " << rain << endl;
  if (1 == rain)
  { // in case of rain read in a file
    param << "Name of the rain file  <rain_NF>:: " << rain_NF << endl;
  }
  else
  {
    param << "Name of the rain file  <rain_NF>:: " << endl;
  }
  param << endl;
  param << "Suffix for the 'Outputs' directory  <suffix_o>:: " << suffix_outputs << endl;
  param << endl;
  if (0 == nbtimes)
  {
    param << "Format of the Output file (1=gnuplot 2=vtk)  <output_f>:: " << endl;
  }
  else
  {
    param << "Format of the Output file (1=gnuplot 2=vtk)  <output_f>:: " << output_format << endl;
  }
  param.close();

  // print the values in the terminal

  cout << "*********************************************************" << endl;
  cout << "number of grid cells Nxcell = " << Nxcell << endl;
  cout << "number of grid cells Nycell = " << Nycell << endl;
  cout << endl;
  cout << "length of the domain L = " << L << endl;
  cout << "length of the domain l = " << l << endl;
  cout << endl;
  cout << "length cells dx = " << dx << endl;
  cout << "length cells dy = " << dy << endl;
  cout << "---------------------------------------------------------" << endl;
  cout << "final time T= " << T << endl;
  cout << "number of times saved= " << nbtimes << endl;

  switch (scheme_type)
  {
  case 1:
    cout << "choice of type of scheme: fixed cfl" << endl;
    cout << "value of the cfl  <cflfix>:: " << cfl_fix << endl;
    break;
  case 2:
    cout << "choice of type of scheme: fixed dt" << endl;
    cout << "timestep (in seconds) <dtfix>:: " << dt_fix << endl;
    cout << "value of the cfl  <cflfix>:: " << cfl_fix << endl;
  }

  cout << "---------------------------------------------------------" << endl;

  switch (Lbound)
  {
  case 1:
    cout << "left condition = imposed h (and q if supercritical) " << endl;
    cout << " imposed discharge  = " << left_imp_discharge << endl;
    cout << " imposed height = " << left_imp_h << endl;
    break;
  case 2:
    cout << "left condition = wall " << endl;
    break;
  case 3:
    cout << "left condition = neumann " << endl;
    break;
  case 4:
    cout << "left condition = periodic " << endl;
    break;
  case 5:
    cout << "left condition = imposed discharge " << endl;
    cout << " imposed discharge = " << left_imp_discharge << endl;
    cout << " imposed height = " << left_imp_h << endl;
  }

  switch (Rbound)
  {
  case 1:
    cout << "right condition = imposed h (and q if supercritical) " << endl;
    cout << " imposed discharge  = " << right_imp_discharge << endl;
    cout << " imposed height = " << right_imp_h << endl;
    break;
  case 2:
    cout << "right condition = wall " << endl;
    break;
  case 3:
    cout << "right condition = neumann " << endl;
    break;
  case 4:
    cout << "right condition = periodic " << endl;
    break;
  case 5:
    cout << "right condition = imposed discharge  " << endl;
    cout << " imposed discharge = " << right_imp_discharge << endl;
    cout << " imposed height = " << right_imp_h << endl;
  }

  switch (Bbound)
  {
  case 1:
    cout << "bottom condition = imposed h (and q if supercritical) " << endl;
    cout << " imposed discharge (inflow supercritical) = " << bottom_imp_discharge << endl;
    cout << " imposed height = " << bottom_imp_h << endl;
    break;
  case 2:
    cout << "bottom condition = wall " << endl;
    break;
  case 3:
    cout << "bottom condition = neumann " << endl;
    break;
  case 4:
    cout << "bottom condition = periodic " << endl;
    break;
  case 5:
    cout << "bottom condition = imposed  discharge " << endl;
    cout << " imposed discharge = " << bottom_imp_discharge << endl;
    cout << " imposed height = " << bottom_imp_h << endl;
  }

  switch (Tbound)
  {
  case 1:
    cout << "top condition = imposed h (and q if supercritical) " << endl;
    cout << " imposed discharge (inflow supercritical) = " << top_imp_discharge << endl;
    cout << " imposed height = " << top_imp_h << endl;
    break;
  case 2:
    cout << "top condition = wall " << endl;
    break;
  case 3:
    cout << "top condition = neumann " << endl;
    break;
  case 4:
    cout << "top  condition = periodic " << endl;
    break;
  case 5:
    cout << "top condition = imposed discharge " << endl;
    cout << " imposed discharge = " << top_imp_discharge << endl;
    cout << " imposed height = " << top_imp_h << endl;
  }

  cout << "---------------------------------------------------------" << endl;

  switch (flux)
  {
  case 1:
    cout << "flux condition = Rusanov" << endl;
    break;
  case 2:
    cout << "flux condition = HLL " << endl;
    break;
  case 3:
    cout << "flux condition = HLL2" << endl;
    break;
  case 4:
    cout << "flux condition = HLLC" << endl;
    break;
  case 5:
    cout << "flux condition = HLLC2" << endl;
  }

  cout << "---------------------------------------------------------" << endl;
  cout << "order of the scheme = " << order << endl;
  cout << "---------------------------------------------------------" << endl;

  switch (fric)
  {
  case 0:
    cout << "no friction" << endl;
    break;
  case 1:
    cout << "friction condition = Manning" << endl;
    if (fric_init == 1)
    {
      cout << "parameter friction = file reading" << endl;
    }
    else
    {
      cout << "parameter friction = " << friccoef << endl;
    }
    break;
  case 2:
    cout << "friction condition = Darcy-Weisbach" << endl;
    if (fric_init == 1)
    {
      cout << "parameter friction = file reading" << endl;
    }
    else
    {
      cout << "parameter friction = " << friccoef << endl;
    }
    break;
  case 3:
    cout << "friction law = laminar" << endl;
    if (fric_init == 1)
    {
      cout << "parameter friction = file reading" << endl;
    }
    else
    {
      cout << "parameter friction = " << friccoef << endl;
    }
  }
  cout << "---------------------------------------------------------" << endl;

  if (2 == order)
  { // in case of order 2
    switch (rec)
    {
    case 1:
      cout << "reconstruction condition = MUSCL" << endl;
      break;
    case 2:
      cout << "reconstruction condition = ENO " << endl;
      cout << "parameter amortENO = " << amortENO << endl;
      break;
    case 3:
      cout << "reconstruction condition = ENO_mod " << endl;
      cout << "parameter amortENO = " << amortENO << endl;
      cout << "parameter modifENO = " << modifENO << endl;
    }

    switch (lim)
    {
    case 1:
      cout << "limitation condition = Minmod " << endl;
      break;
    case 2:
      cout << "limitation condition = VanAlbada " << endl;
      break;
    case 3:
      cout << "limitation condition = VanLeer " << endl;
    }

    cout << "---------------------------------------------------------" << endl;
  }

  switch (inf)
  {
  case 0:
    cout << "no infiltration" << endl;
    break;
  case 1:
    cout << "infiltration model = Green-Ampt" << endl;
    if (Kc_init == 1)
    {
      cout << "parameter Kc     = file reading " << endl;
    }
    else
    {
      cout << "parameter Kc     = " << Kc_coef << endl;
    }
    if (Ks_init == 1)
    {
      cout << "parameter Ks     = file reading " << endl;
    }
    else
    {
      cout << "parameter Ks     = " << Ks_coef << endl;
    }
    if (dtheta_init == 1)
    {
      cout << "parameter dtheta = file reading " << endl;
    }
    else
    {
      cout << "parameter dtheta = " << dtheta_coef << endl;
    }
    if (Psi_init == 1)
    {
      cout << "parameter Psi    = file reading " << endl;
    }
    else
    {
      cout << "parameter Psi    = " << Psi_coef << endl;
    }
    if (zcrust_init == 1)
    {
      cout << "parameter zcrust = file reading " << endl;
    }
    else
    {
      cout << "parameter zcrust = " << zcrust_coef << endl;
    }
    if (imax_init == 1)
    {
      cout << "parameter imax   = file reading " << endl;
    }
    else
    {
      cout << "parameter imax   = " << imax_coef << endl;
    }

    break;
  }

  cout << "---------------------------------------------------------" << endl;
  switch (huv_init)
  {
  case 1:
    cout << "huv initial condition = file reading" << endl;
    break;
  case 2:
    cout << "huv initial condition = h, u & v = 0 " << endl;
    break;
  case 3:
    cout << "huv initial condition = function thacker" << endl;
    break;
  case 4:
    cout << "huv initial condition = function Radial_Dam_dry" << endl;
    break;
  case 5:
    cout << "huv initial condition = function Radial_Dam_wet" << endl;
  }

  cout << "---------------------------------------------------------" << endl;
  switch (rif_init)
  {
    case 1:
      cout << "rif initial condition = file reading" << endl;
      break;
    case 2:
      cout << "rif initial condition = r, i & f = 0 " << endl;
      break;
    case 3:
      cout << "rif initial condition = function thacker" << endl;
      break;
    case 4:
      cout << "rif initial condition = function Radial_Dam_dry" << endl;
    case 5:
      cout << "rif initial condition = function Radial_Dam_wet" << endl;
  }
  
  cout << "---------------------------------------------------------" << endl;
  switch (topo)
  {
  case 1:
    cout << "topography initial condition = file reading" << endl;
    break;
  case 2:
    cout << "topography initial condition =  function generated" << endl;
    break;
  case 3:
    cout << "topography initial condition =  function thacker" << endl;
    break;
  case 4:
    cout << "topography initial condition =  function generated" << endl;
    break;
  case 5:
    cout << "topography initial condition =  function generated" << endl;
  }
  cout << "---------------------------------------------------------" << endl;
  switch (rain)
  {
  case 0:
    cout << "no rain" << endl;
    break;
  case 1:
    cout << "rain condition = file reading" << endl;
    break;
  case 2:
    cout << "rain condition = rain intensity = 0.00001 m/s " << endl;
  }

  switch (output_format)
  {
  case 1:
    cout << "---------------------------------------------------------" << endl;
    cout << "format of the output files = gnuplot " << endl;
    break;
  case 2:
    cout << "---------------------------------------------------------" << endl;
    cout << "format of the output files = vtk " << endl;
    break;
  }
  cout << "---------------------------------------------------------" << endl;
  cout << "entries ok" << endl;
  cout << "*********************************************************" << endl;
}

int Parameters::get_Nxcell() const
{

  /**
   * @details
   * @return The number of cells in space in the first (x) direction Parameters#Nxcell.
   */

  return Nxcell;
}

int Parameters::get_Nycell() const
{

  /**
   * @details
   * @return The number of cells in space in the second (y) direction Parameters#Nycell.
   */

  return Nycell;
}

SCALAR Parameters::get_T() const
{

  /**
   * @details
   * @return The final time Parameters#T.
   */

  return T;
}

int Parameters::get_nbtimes() const
{

  /**
   * @details
   * @return The number of times saved Parameters#nbtimes.
   */

  return nbtimes;
}

int Parameters::get_scheme_type() const
{

  /**
   * @details
   * @return The type of scheme Parameters#scheme_type.
   */

  return scheme_type;
}

SCALAR Parameters::get_dtfix() const
{

  /**
   * @details
   * @return The fixed space step Parameters#dx_fix.
   */

  return dt_fix;
}

SCALAR Parameters::get_cflfix() const
{

  /**
   * @details
   * @return The fixed cfl Parameters#cfl_fix.
   */

  return cfl_fix;
}

SCALAR Parameters::get_dx() const
{

  /**
   * @details
   * @return The space step in the first (x) direction Parameters#dx.
   */

  return dx;
}

SCALAR Parameters::get_dy() const
{

  /**
   * @details
   * @return The space step in the second (y) direction Parameters#dy.
   */

  return dy;
}

int Parameters::get_Lbound() const
{

  /**
   * @details
   * @return The value corresponding to the left boundary condition Parameters#Lbound.
   */

  return Lbound;
}

SCALAR Parameters::get_left_imp_discharge() const
{

  /**
   * @details
   * @return The value of the imposed discharge per cell in the left boundary condition, that is Parameters#left_imp_discharge / Parameters#l.
   */

  return (left_imp_discharge / l);
}

SCALAR Parameters::get_left_imp_h() const
{

  /**
   * @details
   * @return The value of the imposed water height in the left boundary condition Parameters#left_imp_h.
   */

  return left_imp_h;
}

int Parameters::get_Rbound() const
{

  /**
   * @details
   * @return The value corresponding to the right boundary condition Parameters#Rbound.
   */

  return Rbound;
}

SCALAR Parameters::get_right_imp_discharge() const
{

  /**
   * @details
   * @return The value of the imposed discharge per cell in the right boundary condition, that is Parameters#right_imp_discharge / Parameters#l.
   */

  return (right_imp_discharge / l);
}

SCALAR Parameters::get_right_imp_h() const
{

  /**
   * @details
   * @return The value of the imposed water height in the right boundary condition Parameters#right_imp_h.
   */

  return right_imp_h;
}

int Parameters::get_Bbound() const
{

  /**
   * @details
   * @return The value corresponding to the bottom boundary condition Parameters#Bbound.
   */

  return Bbound;
}

SCALAR Parameters::get_bottom_imp_discharge() const
{

  /**
   * @details
   * @return The value of the imposed discharge per cell in the bottom boundary condition, that is Parameters#bottom_imp_discharge / Parameters#L.
   */

  return (bottom_imp_discharge / L);
}

SCALAR Parameters::get_bottom_imp_h() const
{

  /**
   * @details
   * @return The value of the imposed water height in the bottom boundary condition Parameters#bottom_imp_h.
   */

  return bottom_imp_h;
}

int Parameters::get_Tbound() const
{

  /**
   * @details
   * @return The value corresponding to the top boundary condition Parameters#Tbound.
   */

  return Tbound;
}

SCALAR Parameters::get_top_imp_discharge() const
{

  /**
   * @details
   * @return The value of the imposed discharge per cell in the top boundary condition, that is Parameters#top_imp_discharge / Parameters#L.
   */

  return (top_imp_discharge / L);
}

SCALAR Parameters::get_top_imp_h() const
{

  /**
   * @details
   * @return The value of the imposed water height in the bottom boundary condition Parameters#top_imp_h.
   */

  return top_imp_h;
}

int Parameters::get_flux() const
{

  /**
   * @details
   * @return The value corresponding to the flux Parameters#flux.
   */

  return flux;
}

int Parameters::get_order() const
{

  /**
   * @details
   * @return The order of the scheme Parameters#order.
   */

  return order;
}

int Parameters::get_rec() const
{

  /**
   * @details
   * @return The value corresponding to the reconstruction Parameters#rec.
   */

  return rec;
}

int Parameters::get_fric() const
{

  /**
   * @details
   * @return The value corresponding to the friction law Parameters#fric.
   */

  return fric;
}

int Parameters::get_lim() const
{

  /**
   * @details
   * @return The value corresponding to the limiter Parameters#lim.
   */

  return lim;
}

SCALAR Parameters::get_amortENO() const
{

  /**
   * @details
   * @return The value of the amortENO parameter Parameters#amortENO.
   */

  return amortENO;
}

SCALAR Parameters::get_modifENO() const
{

  /**
   * @details
   * @return The value of the modifENO parameter Parameters#modifENO.
   */

  return modifENO;
}

int Parameters::get_inf() const
{

  /**
   * @details
   * @return The value corresponding to the infiltration Parameters#inf.
   */

  return inf;
}

int Parameters::get_fric_init() const
{

  /**
   * @details
   * @return The value corresponding to the friction coefficient Parameters#fric_init.
   */

  return fric_init;
}

string Parameters::get_frictionNameFile(void) const
{

  /**
   * @details
   * @return The friction coefficient path + Input directory Parameters#fric_namefile.
   */

  return fric_namefile;
}

string Parameters::get_frictionNameFileS(void) const
{

  /**
   * @details
   * @return The friction coefficient namefile (inside the Input directory) Parameters#fric_NF.
   */

  return fric_NF;
}

SCALAR Parameters::get_friccoef() const
{

  /**
   * @details
   * @return The value of the friction coefficient Parameters#friccoef.
   */

  return friccoef;
}

string Parameters::get_topographyNameFile(void) const
{

  /**
   * @details
   * @return The topography path + Input directory Parameters#topography_namefile.
   */

  return topography_namefile;
}

string Parameters::get_topographyNameFileS(void) const
{

  /**
   * @details
   * @return The topography namefile (inside the Input directory) Parameters#topo_NF.
   */

  return topo_NF;
}

int Parameters::get_topo() const
{

  /**
   * @details
   * @return The value corresponding to the topography Parameters#topo.
   */

  return topo;
}

int Parameters::get_huv() const
{

  /**
   * @details
   * @return The value corresponding to the initialization of h and u,v Parameters#huv_init.
   */

  return huv_init;
}

int Parameters::get_rif() const
{
  /**
  * @details
  * @return the value corresponding to the initialization of r, i, and f Parameters#rif_init.
  */

  return rif_init;
}

string Parameters::get_huvNameFile(void) const
{

  /**
   * @details
   * @return The h and u,v path for the initialization + Input directory Parameters#huv_namefile.
   */

  return huv_namefile;
}

string Parameters::get_rifNameFile(void) const
{
  /**
  * @details
  * @return the r, i, and f path for the initialization + Input directory Parameters#rif_namefile.
  */

  return rif_namefile;
}

string Parameters::get_huvNameFileS(void) const
{

  /**
   * @details
   * @return The h and u namefile for the initialization (inside the Input directory) Parameters#huv_NF.
   */

  return huv_NF;
}

string Parameters::get_rifNameFileS(void) const
{
  /**
  * @details
  * @return the r, i, and f namefile for the initialization (inside the Input directory) Parameters#rif_NF.
  */
  return rif_NF;

}

int Parameters::get_rain() const
{

  /**
   * @details
   * @return The value corresponding to the initialization of the rain Parameters#rain.
   */

  return rain;
}

string Parameters::get_rainNameFile(void) const
{

  /**
   * @details
   * @return The rain path for the initialization + Input directory Parameters#rain_namefile.
   */

  return rain_namefile;
}

string Parameters::get_rainNameFileS(void) const
{

  /**
   * @details
   * @return The rain namefile for the initialization (inside the Input directory) Parameters#rain_NF.
   */

  return rain_NF;
}

string Parameters::get_outputDirectory(void) const
{

  /**
   * @details
   * @return The output directory with the suffix Parameters#output_directory.
   */

  return output_directory;
}

int Parameters::get_Kc_init() const
{

  /**
   * @details
   * @return The value corresponding to the initialization of Kc Parameters#Kc_init.
   */

  return Kc_init;
}

SCALAR Parameters::get_Kc_coef() const
{

  /**
   * @details
   * @return The value of Kc Parameters#Kc_coef.
   */

  return Kc_coef;
}

string Parameters::get_KcNameFile(void) const
{

  /**
   * @details
   * @return The Kc path for the initialization + Input directory Parameters#Kc_namefile.
   */

  return Kc_namefile;
}

string Parameters::get_KcNameFileS() const
{

  /**
   * @details
   * @return The Kc namefile for the initialization (inside the Input directory) Parameters#Kc_NF.
   */

  return Kc_NF;
}

int Parameters::get_Ks_init() const
{

  /**
   * @details
   * @return The value corresponding to the initialization of Ks Parameters#Ks_init.
   */

  return Ks_init;
}

SCALAR Parameters::get_Ks_coef() const
{

  /**
   * @details
   * @return The value of Ks Parameters#Ks_coef.
   */

  return Ks_coef;
}

string Parameters::get_KsNameFile(void) const
{

  /**
   * @details
   * @return The Ks path for the initialization + Input directory Parameters#Ks_namefile.
   */

  return Ks_namefile;
}

string Parameters::get_KsNameFileS() const
{

  /**
   * @details
   * @return The Ks namefile for the initialization (inside the Input directory) Parameters#Ks_NF.
   */

  return Ks_NF;
}

int Parameters::get_dtheta_init() const
{

  /**
   * @details
   * @return The value corresponding to the initialization of dtheta Parameters#dtheta_init.
   */

  return dtheta_init;
}

SCALAR Parameters::get_dtheta_coef() const
{

  /**
   * @details
   * @return The value of dtheta Parameters#dtheta_coef.
   */

  return dtheta_coef;
}

string Parameters::get_dthetaNameFile(void) const
{

  /**
   * @details
   * @return The dtheta path for the initialization + Input directory Parameters#dtheta_namefile.
   */

  return dtheta_namefile;
}

string Parameters::get_dthetaNameFileS() const
{

  /**
   * @details
   * @return The dtheta namefile for the initialization (inside the Input directory) Parameters#dtheta_NF.
   */

  return dtheta_NF;
}

int Parameters::get_Psi_init() const
{

  /**
   * @details
   * @return The value corresponding to the initialization of Psi Parameters#Psi_init.
   */

  return Psi_init;
}

SCALAR Parameters::get_Psi_coef() const
{

  /**
   * @details
   * @return The value of Psi Parameters#Psi_coef.
   */

  return Psi_coef;
}

string Parameters::get_PsiNameFile(void) const
{

  /**
   * @details
   * @return The Psi path for the initialization + Input directory Parameters#Psi_namefile.
   */

  return Psi_namefile;
}

string Parameters::get_PsiNameFileS() const
{

  /**
   * @details
   * @return The Psi namefile for the initialization (inside the Input directory) Parameters#Psi_NF.
   */

  return Psi_NF;
}

int Parameters::get_zcrust_init() const
{

  /**
   * @details
   * @return The value corresponding to the initialization of zcrust Parameters#zcrust_init.
   */

  return zcrust_init;
}

SCALAR Parameters::get_zcrust_coef() const
{

  /**
   * @details
   * @return The value of zcrust Parameters#zcrust_coef.
   */

  return zcrust_coef;
}

string Parameters::get_zcrustNameFile(void) const
{

  /**
   * @details
   * @return The zcrust path for the initialization + Input directory Parameters#zcrust_namefile.
   */

  return zcrust_namefile;
}

string Parameters::get_zcrustNameFileS() const
{

  /**
   * @details
   * @return The zcrust namefile for the initialization (inside the Input directory) Parameters#zcrust_NF.
   */

  return zcrust_NF;
}

int Parameters::get_imax_init() const
{

  /**
   * @details
   * @return The value corresponding to the initialization of imax Parameters#imax_init.
   */

  return imax_init;
}

SCALAR Parameters::get_imax_coef() const
{

  /**
   * @details
   * @return The value of imax Parameters#imax_coef.
   */

  return imax_coef;
}

string Parameters::get_imaxNameFile(void) const
{

  /**
   * @details
   * @return The imax path for the initialization + Input directory Parameters#imax_namefile.
   */

  return imax_namefile;
}

string Parameters::get_imaxNameFileS() const
{

  /**
   * @details
   * @return The imax namefile for the initialization (inside the Input directory) Parameters#imax_NF.
   */

  return imax_NF;
}

string Parameters::get_suffix(void) const
{

  /**
   * @details
   * @return The suffix (for the output directory) Parameters#suffix_outputs.
   */

  return suffix_outputs;
}

int Parameters::get_output() const
{

  /**
   * @details
   * @return The type of output Parameters#output_format.
   */

  return output_format;
}

void Parameters::fill_array(TAB &myarray, const SCALAR myvalue) const
{

  /**
   * @details Fills an array with a constant value.
   * @param[in, out] myarray array to fill.
   * @param[in] myvalue value.
   */

  for (unsigned int i = 0; i < myarray.size(); i++)
  {
    //fill is a function provided by the STL that initializes a container with a value.
    fill(myarray[i].begin(), myarray[i].end(), myvalue);
  } //end for i
}

void Parameters::fill_array(TAB &myarray, string namefile) const
{

  /**
   * @details Fills an array with the values given in the file
   * @param[in, out] myarray array to fill.
   * @param[in] namefile name of the file containing the values to be inserted into the array.
   * @warning ***: ERROR: cannot open the file.
   * @warning ***: ERROR: the number of data in this file is too big/small.
   * @warning ***: ERROR: line ***.
   * @warning ***: ERROR: the value for the point x = *** y = *** is missing.
   * @warning ***: WARNING: line *** ; a commentary should begin with the # symbol.
   * @note If the array cannot be filled correctly, the code will exit with failure termination code.
   */

  SCALAR x, y, value;
  int row, column;
  int i = 0;       //iterator
  int num_lin = 0; //iterator
  string line;     // string to store a line of the input file
  char car;        //First character of a commentary line

  //Opening the file and verification if it exists.
  ifstream getdata(namefile.c_str(), ios::in);   //ios::in open for input operations
  if (!getdata)
  {
    //The input file cannot be opened.
    cerr << namefile << ": ERROR: cannot open the file\n";
    exit(EXIT_FAILURE);
  }

  //Initialization of myarray to the largest finite representable floating-point number.
  //So that we will be able to check that the user fills the variable properly.
  for (int i = 1; i < Nxcell + 1; i++)
  {
    for (int j = 1; j < Nycell + 1; j++)
    {
      myarray[i][j] = MAX_SCAL;
    } //end for j
  }   //end for i

  // As long as the end of the file is not reached, the next line is read.
  while (!getdata.eof())
  {
    num_lin++;
    getline(getdata, line, '\n'); // read a line  istream& getline (istream& is,string& str, char delim)
    //Extract characters from "is" and stores them into str until the delimitation character delim is found
    istringstream entree(line);
    if (entree >> x >> y >> value)
    {

      //Check that the input file contains the expected number of data.
      if (++i > Nxcell * Nycell)
      {
        cerr << namefile << ": ERROR: the number of data in this file is too big!" << endl;
        exit(EXIT_FAILURE);
      }
      //We compute the index of the array from the space variable in the input file.
      row = (int)ceil(x / dx);    //Rounds x upward, returing the smallest integral value that is not less than x
      column = (int)ceil(y / dy); //The index corresponds to the smallest integer superior or equal to x/dx or y/dy.

      //Error if the x or y value of the input file is out of the domain.
      if (x < 0)
      {
        cerr << namefile << ": ERROR: line " << num_lin << ", x = " << x << " must be positive." << endl;
        exit(EXIT_FAILURE);
      }
      if (y < 0)
      {
        cerr << namefile << ": ERROR: line " << num_lin << ", y = " << y << " must be positive." << endl;
        exit(EXIT_FAILURE);
      }
      if (x > Nxcell * dx)
      {
        cerr << namefile << ": ERROR: line " << num_lin << ", x = " << x << " must be lower than " << Nxcell * dx << "." << endl;
        exit(EXIT_FAILURE);
      }
      if (y > Nycell * dy)
      {
        cerr << namefile << ": ERROR: at line " << num_lin << ", y = " << y << " must be lower than " << Nycell * dy << "." << endl;
        exit(EXIT_FAILURE);
      }

      //Error if the x or y value of the input file is not near the center of the cell.
      if (fabs(x - (row - 0.5) * dx) > dx * RATIO_CLOSE_CELL)
      {
        cout << namefile << ": ERROR: line " << num_lin << "; x = " << x << ";\n This value is not close enough to the center of the cell.\n You may want to replace it with " << (row - 0.5) * dx << endl;
        exit(EXIT_FAILURE);
      }
      if (fabs(y - (column - 0.5) * dy) > dy * RATIO_CLOSE_CELL)
      {
        cout << "column =" << column << endl;
        cout << namefile << ": ERROR: line " << num_lin << "; y = " << y << ";\n This value is not close enough to the center of the cell.\n You may want to replace it with " << (column - 0.5) * dy << endl;
        exit(EXIT_FAILURE);
      }

      //Store the input values into the array.
      myarray[row][column] = value;
    }
    else
    {
      car = '#'; //Initialization of the character used to identity the beginning of a comment
      istringstream entree_car(line);
      entree_car >> car;
      if (car != '#')
      {
        cout << namefile << ": WARNING: line " << num_lin << "; a commentary should begin with the # symbol " << endl;
      }
    }
  }

  //Closing the input file
  getdata.close();

  //Check that the input file contains the expected number of data.
  if (i < Nxcell * Nycell)
  {
    cerr << namefile << ": ERROR: the number of data in this file is too small!" << endl;
    exit(EXIT_FAILURE);
  }

  //Final check: Does all the grid cells were filled with a value?
  for (int i = 1; i < Nxcell + 1; i++)
  {
    for (int j = 1; j < Nycell + 1; j++)
    {
      if (myarray[i][j] >= MAX_SCAL)
      {
        cerr << namefile << ": ERROR: the value for the point x =" << (i - 0.5) * dx << " y = " << (j - 0.5) * dy << " is missing!" << endl;
        exit(EXIT_FAILURE);
      }
    } //end for j
  }   //end for i
}

Parameters::Parameters()
{
} //construntor

Parameters::~Parameters()
{
} //destructor
