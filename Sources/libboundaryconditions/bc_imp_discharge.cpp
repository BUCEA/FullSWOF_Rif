/**
 * @file bc_imp_discharge.cpp
 * @author Ulrich Razafison <ulrich.razafison@math.cnrs.fr> (2011)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @version 1.06.01
 * @date 2015-10-29
 *
 * @brief Imposed discharge
 * @details
 * Boundary condition:
 * imposed discharge (and water height if necessary).
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

#include "bc_imp_discharge.hpp"

Bc_imp_discharge::Bc_imp_discharge(Parameters & par,TAB & z,int n1, int n2): Boundary_condition(par){

  /**
   * @details
   * @param[in] par parameter, contains all the values from the parameters file (unused).
   * @param[in] n1 integer to specify whether it is the left (-1) or the right (1) boundary.
   * @param[in] n2 integer to specify whether it is the bottom (-1) or the top (1) boundary.
   * @param[in,out] z vector that represents the topography with suitable values on the fictive cells.
   */
  
  flag=0;
  
  tol = 0.1*HE_CA;
  maxiter = 10000;
  unormbound = 0.;
  utanbound = 0.;
  hbound = 0.;
  
  // Initialization of the topography at the boundary
  
  if ((-1==n1)&&(0==n2)) {
    // Left Boundary condition
    for (int j=1 ; j<NYCELL+1 ; j++){
      z[0][j]=2*z[1][j]-z[2][j];
    }
  }
  if ((0==n1)&&(-1==n2)) {
    // Bottom Boundary condition
    for (int i=1 ; i<NXCELL+1 ; i++){
      z[i][0]=2*z[i][1]-z[i][2];
    }
  }
  
  if ((1==n1)&&(0==n2)) {
    //Right  Boundary condition
    for (int j=1 ; j<=NYCELL ; j++){
      z[NXCELL+1][j]=2*z[NXCELL][j]-z[NXCELL-1][j];
    }
  }
  
  if ((0==n1)&&(1==n2)) {
    // Top Boundary condition
    for (int i=1 ; i<NXCELL+1 ; i++){
      z[i][NYCELL+1]=2*z[i][NYCELL]-z[i][NYCELL-1];
    }
  }
}

SCALAR Bc_imp_discharge::getValueOfPolynomial(const SCALAR HIN, const SCALAR UNORM_IN, const SCALAR UTAN_IN, const SCALAR QFIX, const SCALAR H,int n1, int n2) const{
  
  /**
   * @details
   * Computes @f$2H\sqrt{gH}-(n1+n2)(UNORM\_IN+2(n1+n2)\sqrt{gHIN})H-|QFIX|@f$ where @f$n1, n2@f$ are the normals.
   * @param[in] HIN water height of the first cell inside the domain.
   * @param[in] UNORM_IN normal velocity of the first cell inside the domain.
   * @param[in] UTAN_IN tangential velocity of the first cell inside the domain (unused).
   * @param[in] QFIX fixed (imposed) value of the discharge.
   * @param[in] H value for the variable of the polynomial function.
   * @param[in] n1 integer to specify whether it is the left (-1) or the right (1) boundary.
   * @param[in] n2 integer to specify whether it is the bottom (-1) or the top (1) boundary.
   * @return The value of the polynomial function.
   */
  
  (void) UTAN_IN; //unused variable

  return 2 * sqrt(GRAV*H) * H -(n1+n2)*(UNORM_IN +(n1+n2)*2*sqrt(GRAV * HIN)) * H -fabs(QFIX);//fabs compute absolute value of QFIX
}

SCALAR Bc_imp_discharge::getValueofDerivativeOfPolynomial(const SCALAR HIN, const SCALAR UNORM_IN, const SCALAR UTAN_IN, const SCALAR H,int n1, int n2) const{
  
  /**
   * @details
   * Computes @f$ 3\sqrt{gH}-(n1+n2)(UNORM\_IN+2(n1+n2)\sqrt{gHIN})@f$ where @f$n1, n2@f$ are the normals.
   * @param[in] HIN water height of the first cell inside the domain.
   * @param[in] UNORM_IN normal velocity of the first cell inside the domain.
   * @param[in] UTAN_IN tangential velocity of the first cell inside the domain (unused).
   * @param[in] H value for the variable of the polynomial function.
   * @param[in] n1 integer to specify whether it is the left (-1) or the right (1) boundary.
   * @param[in] n2 integer to specify whether it is the bottom (-1) or the top (1) boundary.
   * @return The value of derivative of the polynomial function defined in Bc_imp_discharge::getValueOfPolynomial().
   */

  (void) UTAN_IN; //unused variable
  
  return 3 * sqrt(GRAV*H)-(n1+n2)*(UNORM_IN +(n1+n2)*2*sqrt(GRAV * HIN));
}

SCALAR Bc_imp_discharge::newtonSolver(const SCALAR HIN, const SCALAR UNORM_IN, const SCALAR UTAN_IN, const SCALAR QFIX, const SCALAR H_INIT,int n1, int n2) const{
  
  /**
   * @details
   * Finds the root of the polynomial function corresponding to the imposed discharge.
   * Needs the evaluation of the function and of its derivative.
   * @param[in] HIN water height of the first cell inside the domain.
   * @param[in] UNORM_IN normal velocity of the first cell inside the domain.
   * @param[in] UTAN_IN tangential velocity of the first cell inside the domain.
   * @param[in] QFIX fixed (imposed) value of the discharge.
   * @param[in] H_INIT initialization of the Newton solver.
   * @param[in] n1 integer to specify whether it is the left (-1) or the right (1) boundary.
   * @param[in] n2 integer to specify whether it is the bottom (-1) or the top (1) boundary.
   * @warning Warning: Newton bc did not converge.
   * @return <em> \b h</em>: water height that satifies Riemann invariants.
   */
  
  SCALAR norm, F, gradF, h,sigma;
  int iter = 0;
  // initialization
  h = H_INIT;
  F = getValueOfPolynomial(HIN,UNORM_IN,UTAN_IN,QFIX,h,n1,n2);
  norm = fabs(F);
  while ((norm > tol) && (iter <= maxiter)){
    gradF = getValueofDerivativeOfPolynomial(HIN,UNORM_IN,UTAN_IN,h,n1,n2);
    sigma = -F / gradF;
    h += sigma;
    F = getValueOfPolynomial(HIN,UNORM_IN,UTAN_IN,QFIX,h,n1,n2);
    norm = fabs(F);
    iter += 1;
  }// end while
  if (maxiter+1 == iter) {
    cout << "Warning: Newton bc did not converge" << endl;
  }
  
  return h;
  
}

void Bc_imp_discharge::calcul(SCALAR hin, SCALAR unorm_in, SCALAR utan_in, SCALAR hfix, SCALAR qfix, SCALAR hin_oppbound, SCALAR unorm_in_oppbound, SCALAR utan_in_oppbound, SCALAR time,int n1, int n2){
  
  /**
   * @details
   * Two cases are considered: subcritical and supercritical flows.
   * @param[in] hin water height of the first cell inside the domain.
   * @param[in] unorm_in normal velocity of the first cell inside the domain.
   * @param[in] utan_in tangential velocity of the first cell inside the domain.
   * @param[in] hfix fixed (imposed) value of the water height (only for the supercritical case).
   * @param[in] qfix fixed (imposed) value of the discharge.
   * @param[in] hin_oppbound value of the water height of the first cell inside the domain at the opposite bound (unused).
   * @param[in] unorm_in_oppbound value of the normal velocity of the first cell inside the domain at the opposite bound (unused).
   * @param[in] utan_in_oppbound value of the tangential velocity of the first cell inside the domain at the opposite bound (unused).
   * @param[in] time current time (unused).
   * @param[in] n1 integer to specify whether it is the left (-1) or the right (1) boundary.
   * @param[in] n2 integer to specify whether it is the bottom (-1) or the top (1) boundary.
   * @warning Warning in the method Bc_imp_discharge::calcul() The water height at the inflow is zero ... continuing!
   * @par Modifies
   * Boundary_condition#hbound water height on the fictive cell.\n
   * Boundary_condition#unormbound normal velocity on the fictive cell.\n
   * Boundary_condition#utanbound tangential velocity on the fictive cell.
   */
  
  (void) hin_oppbound; //unused variable
  (void) unorm_in_oppbound; //unused variable
  (void) utan_in_oppbound; //unused variable
  (void) time; //unused variable
  
  
  if (fabs(unorm_in)-sqrt(GRAV*hin)<=0.) { // sub critical (fluvial)
    h_init = 2.0*pow(fabs(qfix)/sqrt(GRAV), (SCALAR) 2./3.)+1.0;
    // we initialize the water height with h_init > h_critical = (|qfix|/sqrt(g))^2/3
    // in order to converge to the root we want
    h = newtonSolver(hin,unorm_in,utan_in,qfix,h_init,n1,n2);
    if (h <= HE_CA){
      if (0 == flag ){
        cerr << "Warning in the method bc_imp_discharge::calcul()\n";
        cerr << "The water height at the inflow is zero ... continuing! \n";
        flag =1;
      }
      unormbound = 0.;
      utanbound = 0.;
      hbound = 0.;
    }// end if
    else{
      unormbound = qfix / h;
      utanbound = 0;
      hbound = h;
    }// end else
  
  }else { // super critical (torrential)
    if (hfix < HE_CA) {
      unormfix = 0.;
      hfix = 0.;
    }
    else {
      unormfix = qfix/hfix;
    }
    hbound = hfix;
    unormbound = unormfix;
    utanbound = 0.;
  }
}


Bc_imp_discharge::~Bc_imp_discharge(){
}
