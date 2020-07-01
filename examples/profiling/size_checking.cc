/**
 * @file   size_checking.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Feb 22 15:28:23 2018
 *
 * @brief  Script used to predict the memory required to run a simulation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "spectral_model.hh"
#include "simulation_driver.hh"
#include "interfacer.hh"
#include "cohesive_law.hh"
#include "coulomb_law.hh"
#include "regularized_coulomb_law.hh"
#include "data_dumper.hh"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <string>
#include <iomanip>
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]){

  // Note : Construct the pre-integrated material kernels before running this simulation
  // Use "invert_serial.f" to construct kernel files

  std::cout << "./size_checking <nb_ele_x> <nb_ele_z> [ <nb_time_steps>=0 <beta=0.2> <dom_sizex>=1.0 <dom_sizez>=1.0]" << std::endl;

  UInt nex = std::atoi(argv[1]); 
  UInt nez = std::atoi(argv[2]);    
  UInt nb_time_steps = 0;
  Real beta = 0.2;
  Real dom_sizex = 1.0;
  Real dom_sizez = 1.0;
  
  // Geometry description
  if(argc>3)
    nb_time_steps = std::atoi(argv[3]); 
  if(argc>4)
    beta = std::atof(argv[4]); 
  if(argc>5)
    dom_sizex = std::atof(argv[5]); 
  if(argc>6)
    dom_sizez = std::atof(argv[6]);
   
  Real nu =  0.35;
  Real E = 81e9;
  Real cs = 3000;

  // Cut of the loaded material kernels
  UInt tcut = 100; 
  
  SpectralModel model({nex,nez}, nb_time_steps, {dom_sizex,dom_sizez},
		      nu, nu, E, E, cs, cs, tcut, tcut,
		      "Blank simulation used to estimate the required memory");

  model.initModel(beta,true);
       
  return 0;
}
