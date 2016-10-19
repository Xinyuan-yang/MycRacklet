/**
 * @file   single_mat_heterog.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Nov 11 09:19:07 2015
 *
 * @brief  Homogeneous fracture at a constant crack speed
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

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]){

  // Note : Construct the pre-integrated material kernels before running this simulation
  // Use "invert_serial.f" to construct kernel files

  // Required command line arguments. [...] delimitates the optional ones
  std::cout << " ./cst_speed_homog <crack_speed> <loading_file> <load_writing?> [ <output_folder_name>=./ <nb_t_steps>=10000 <nb_ele>=4096 ]" << std::endl;

  std::string sim_name = "Single material homogeneous interface controled rupture speed";

  // Geometry description
  UInt nb_time_steps = 10000; 
  UInt nb_elements = 4096;
  Real dom_size = 1.0;
  Real crack_size = 0.05;
  Real nu =  0.35;
  Real E = 5.3e9;
  Real cs = 1263;

  // Cut of the loaded material kernels
  UInt tcut = 100;
  
  // Loading case
  Real load = 3e6;
  Real psi = 90;
  Real phi = 0;
  UInt l_index = 1;

  // Target crack speed v=cr_speed*c_s
  Real cr_speed = std::atof(argv[1]);
  // Name of the loading file
  std::string load_file = argv[2];
  // True=simulation to tailor the loading file
  // False=simulation reading the loading file
  bool write = (bool)(std::atoi(argv[3]));

  // Output folder
  std::string output_folder = "./";
  if(argc > 4)
    output_folder = argv[4];

  if(argc > 5)
    nb_time_steps = std::atoi(argv[5]);

  if(argc > 6)
    nb_elements = std::atoi(argv[6]);

  std::cout << "./cst_speed_homog " 
	    << "v_imposed:" << cr_speed << " "
	    << "load_file:" << load_file << " "
	    << "output folder:" << output_folder << " " 
	    << "nb_time steps:" << nb_time_steps << " " 
	    << "nb_elements:" << nb_elements << " "
	    << std::endl;
  
  // Cohesive paramters
  Real crit_n_open = 0.02e-3;
  Real crit_s_open = 0.02e-3;
  Real max_s_str = 5e6;
  Real max_n_str = 5e6;

  FractureLaw * fracturelaw; 

  // Friction paramters
  bool overlap = 0;
  Real regularized_time_scale = 0.1;
  Real coef_frict = 0.25;
  ContactLaw * contactlaw = new RegularizedCoulombLaw(coef_frict, regularized_time_scale, nb_elements);

  /* -------------------------------------------------------------------------- */

  SpectralModel model({nb_elements,1}, nb_time_steps, {dom_size,0.}, 
		      nu, nu, E, E, cs, cs, tcut, tcut, overlap, l_index, 
		      fracturelaw, contactlaw, sim_name, output_folder); 

  // SimulationDriver object helping to launch controlled-speed simulation
  SimulationDriver sim_driver(model, cr_speed, dom_size/2);
  
  Interfacer<_linear_coupled_cohesive> interfacer(model);
  interfacer.createThroughCenteredCrack(crack_size, crit_n_open, crit_s_open, 
					max_n_str, max_s_str);
  interfacer.applyInterfaceCreation();

  if(write){
    // Init algorithm to tailor constant speed loading conditions
    sim_driver.initConstantSpeed(load, psi, phi, max_s_str);
  }
  else {
    // Read an existing file to set loading conditions
    sim_driver.initLoadingFromFile(output_folder+load_file);
  }
    
  DataDumper dumper(model);

  dumper.initEnergetics("Energy.cra");
  dumper.initDumper("ST_Diagram_id.cra", _id_crack, 0.5, 1, 0.5*nb_elements);
  dumper.initDumper("ST_Diagram_shear_velo_jump.cra", _shear_velocity_jumps, 0.5, 1, 0.5*nb_elements);

  UInt x_tip;
   
  for (UInt t = 0; t < nb_time_steps ; ++t) {

    // High level method embedding the resolution of one time step by the SpectralModel
    x_tip = sim_driver.solveStep();  

    if(x_tip>0.95*nb_elements)
      break;

    if(t%5==0){
      dumper.dumpAll();
      dumper.printEnergetics();
    }

    if (t%((UInt)(0.1*nb_time_steps))==0){
      std::cout << "Process at " << (Real)t/(Real)nb_time_steps*100 << "% " 
		<< " crack tip at " << (Real)(x_tip)/(Real)nb_elements*100 << "% "
		<< std::endl;
    }
  }

  // Do not forget to create the resulting loading file
  if (write)
    sim_driver.writeLoading(output_folder+load_file);

  delete contactlaw;
  
  return 0;
}
