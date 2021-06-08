/**
 * @file   test_multimodal_convolution
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Dec 7 17:56:13 2015
 *
 * @brief  Test computing a discrete convolution integral
 *
 * @section LICENSE
 *
 * cRacklet - A spectral boundary integral method for interface fracture simulation
 * Copyright (©) 2012 - 2013 Fabian Barras
 *               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 *               Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 * 
 * cRacklet is the result of a collaboration between the Computational Solid Mechanics 
 * Laboratory (LSMS) of Ecole Polytechnique Fédérale de Lausanne (EPFL), Switzerland 
 * and the Department of Aerospace Engineering of the University of Illinois at 
 * Urbana-Champaign, United States of America.
 * 
 * cRacklet is free software: you can redistribute it and/or modify it under the terms 
 * of the GNU General Public License as published by the Free Software Foundation, 
 * either version 3 of the License, or (at your option) any later version.
 * 
 * cRacklet is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program.  
 * If not, see <http://www.gnu.org/licenses/>.
 */
/* -------------------------------------------------------------------------- */
#include "spectral_convolution_manager.hh"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <gsl/gsl_sf_bessel.h>
#include <assert.h>
#include <math.h>
/* -------------------------------------------------------------------------- */
struct H33 : public ConvolutionManager::KernelFunctor{
  H33(){};
  virtual ~H33(){};
  Real operator()(Real tau){

    if(tau!=0)
      return gsl_sf_bessel_J1(tau)/tau;
    else
      return 0.5;
  }
};

int main() {

  UInt nb_time = 200000;
  UInt nb_el = 8;
  Real X = 1.;
  Real dx = X/(Real)nb_el;
  Real q0 = 2*M_PI/X;
  Real E = 5e9;
  Real nu = 0.35;
  Real mu  = E/(2*(1+nu));
  Real cs = 1200;
  Real cut = 100.;
  Real beta = 0.4;
  
  UInt nb_modes = nb_el/2;
  std::vector<Real> j_ksi(nb_modes);

  for (UInt i=0 ; i < nb_modes; ++i) {
    j_ksi[i] = q0*Real(i+1);
  }

  Real dt = beta*dx;

  std::vector<Real> t(nb_time);

  SpectralConvolutionManager convo_manager(dt, j_ksi, 1, 1);

  convo_manager.init(cut, nb_time,false);

  ConvolutionManager::KernelFunctor ** h33 = convo_manager.getKernelFunctor(0);
  (*h33) = new H33;

  std::cout.precision(6);
  UInt wdth = 15;

  Real * U = new Real[2*nb_modes];
  Real * res = new Real[2*nb_modes];
 
  for(UInt i = 0; i < nb_time; ++i) {

    t[i] = (Real)(dt/cs)*i;
    Real * it = U; 
    for (UInt j=0 ; j < nb_modes; ++j) {
      for (UInt k=0; k<2; ++k){
	*it = exp(-10*pow(t[i]-2.5,2));
	++it;
      }
    }

    std::cout << std::setw(wdth) << *U;
    convo_manager.storeFields(U);
    convo_manager.computeConvolution(res,0,0);

    Real * resit = res;
    for (UInt j=0 ; j < nb_modes; ++j) {
      std::cout << std::setw(wdth) << mu*j_ksi[j]*(*resit);
      resit+=2;
    } 
    std::cout << std::endl; 
  }
  
  delete[] U;
  delete[] res;
  return 0;
}
