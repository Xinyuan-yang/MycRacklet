/**
 * @file   printKernels.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Mar 16 17:07:14 2017
 *
 * @brief  Scipt to output kernel function
 *
 * @section LICENSE
 *
 * cRacklet - A spectral boundary integral method for interface fracture simulation
 * Copyright (©) 2012 - 2013 Fabian Barras
 *               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
int main() {

  UInt nb_time = 50000;
  UInt nb_el = 2;
  Real X = 1.;
  Real dx = X/(Real)nb_el;
  Real q0 = 2*M_PI/X;
  Real E = 5e9;
  Real nu = 0.35;
  Real mu  = E/(2*(1+nu));
  Real cs = 1200;
  Real cut = 100.;
  Real beta = 0.4;
  UInt nt_cut = std::min(nb_time,
			   (UInt)(cut*(Real)nb_el*0.5/(M_PI*beta))-1);
  std::cout << "Nb stored time step " << nt_cut << std::endl;
  Real dt = q0*beta*dx;

  std::vector<Real> t(nb_time);

  ConvolutionManager convo_manager(dt);

  ConvolutionManager::H_33 * h33 = new ConvolutionManager::H_33;
  
  convo_manager.init(nt_cut);
  convo_manager.loadKernel(h33);
  convo_manager.dumpKernel("kernel_33.txt",0.,25.,250);
  std::cout << "Kernel printed !" << std::endl;

  return 0;
}
