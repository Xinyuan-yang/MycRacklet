/**
 * @file   test_constant_speed_writing_space_ctrl.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Feb 15 11:04:01 2018
 *
 * @brief  Testing constant speed fracture (space-controlled loading)
 *
 * @section LICENSE
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
 */
/* -------------------------------------------------------------------------- */
#include "spectral_model.hh"
#include "simulation_driver.hh"
#include "interfacer.hh"
#include "cohesive_law.hh"
#include "data_dumper.hh"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <string>

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]){

  // Note : Construct the pre-integrated material kernels before running this simulation
  // Use "invert_serial.f" to construct kernel files

  std::string sim_name = "Testing controlled rupture speed";

  // Geometry description
  UInt nb_elements = 512;
  Real nu =  0.33;
  Real E = 5.3e9;
  Real cs = 1263;

  // Cut of the loaded material kernels
  UInt tcut = 100;
  
  // Loading case
  Real load = 3e6;
  Real psi = 0;
  Real phi = 0;

  // Target crack speed v=cr_speed*c_s
  Real cr_speed = 0.5;
  // Name of the loading file
  std::string load_file = "test_time_loading_file.txt";

  std::vector<bool> write = {true,false};

  // Output folder
  std::string output_folder = "./";

  // Cohesive paramters
  Real crit_n_open = 0.02e-3;
  Real crit_s_open = 0.02e-3;
  Real max_n_str = 9e6;
  Real max_s_str = 9e6;

  Real G_length = crit_s_open*max_s_str/(load*load*M_PI)*E/(1-nu*nu);

  Real dom_size = 10*G_length;
  Real dx = dom_size/double(nb_elements);
  Real wall_position = 0.95*dom_size;
  UInt propagation_domain = 0.75*nb_elements;
  
  std::cout.precision(6);
  
  std::cout << "Griffith length: " << G_length << std::endl; 
  
  /* -------------------------------------------------------------------------- */

  for (std::vector<bool>::iterator it = write.begin(); it < write.end(); ++it) {

    SpectralModel * model = new SpectralModel (nb_elements, 0, dom_size, 
					       nu, E, cs, tcut,
					       sim_name, output_folder);
    
    // SimulationDriver object helping to launch controlled-speed simulation
    SimulationDriver sim_driver(*model, cr_speed, 0.0);
  
    Interfacer<_linear_coupled_cohesive> interfacer(*model);

    DataRegister::registerParameter("critical_normal_opening",crit_n_open);
    DataRegister::registerParameter("critical_shear_opening",crit_s_open);
    DataRegister::registerParameter("max_normal_strength",max_n_str);
    DataRegister::registerParameter("max_shear_strength",max_s_str);
    interfacer.createUniformInterface();
    
    interfacer.createThroughCrack(0.,5*dx);
    interfacer.createThroughWall(wall_position,dom_size);

    if(*it){
      // Init algorithm to tailor constant speed loading conditions
      sim_driver.initConstantSpeed(load, psi, phi, max_s_str, 0.0, _space_control, 0.9, G_length);
    }
    else {
      // Read an existing file to set loading conditions
      sim_driver.initLoadingFromFile(output_folder+load_file, _space_control, load);
    }
    
    UInt x_tip=0;

    sim_driver.launchCrack(0.,2*G_length,0.05);

    UInt t = 0;

    UInt max_t_step = 1.5*propagation_domain/(model->getParameter("beta")*cr_speed);

    UInt print_info = 0.05*max_t_step;

    const CrackProfile * loading = model->readData(_top_loading);
    UInt wdth = 15;
    
    while ((x_tip<(propagation_domain))&&(t<max_t_step)) {

      // High level method embedding the resolution of one time step by the SpectralModel
      x_tip = sim_driver.solveStep();  

      if (t%((UInt)print_info)==0){
	std::cout << "Crack tip at "
		  << (Real)(x_tip)/(Real)nb_elements*100 << "% " << " and far-field loading "
		  << std::setw(wdth) << (*loading)[0] <<" "
		  << std::setw(wdth) << (*loading)[1] <<" "
		  << std::setw(wdth) << (*loading)[2]
		  << std::endl;
      }
      ++t;
    }
    // Do not forget to create the resulting loading file
    if (*it)
      sim_driver.writeLoading(output_folder+load_file);
    
    delete model;
  }
  return 0;
}
