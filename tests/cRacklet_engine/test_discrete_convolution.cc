/**
 * @file   test_discrete_convolution.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Nov 12 17:56:13 2015
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
#include "convolution_manager.hh"
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

  UInt nb_time = 8192;
  UInt t_cut = 3000;
  Real dt = 10./(Real)nb_time;
  std::vector<Real> t(nb_time);

  H33 * h33 = new H33;

  ConvolutionManager convo_manager(5*dt);
  
  convo_manager.init(t_cut);
  convo_manager.loadKernel(h33);

  std::cout.precision(6);
  UInt wdth = 15;

  Real U;

  for(UInt i = 0; i < nb_time; ++i) {

    t[i] = dt*i;
    U = exp(-10*pow(t[i]-2.5,2));
    convo_manager.storeFields(U);
     
    std::cout << std::setw(wdth) << U 
	      << std::setw(wdth) << convo_manager.computeConvolution() 
	      << std::endl; 
  }
  
  return 0;
}
