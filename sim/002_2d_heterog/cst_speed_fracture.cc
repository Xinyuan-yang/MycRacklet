/**
 * @file   cst_speed_fracture.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Nov 11 09:19:07 2015
 *
 * @brief  Dynamic fracture at a constant crack speed in presence of heterogeneities
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
#include <memory>
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]){

  // Note : Construct the pre-integrated material kernels before running this simulation
  // Use "invert_serial.f" to construct kernel files
  // Required command line arguments. [...] delimitates the optional ones

  std::cout << " ./cst_speed_fracture <crack_speed> <loading_file> <load_writing?> [ <output_folder_name>=./ <nb_t_steps>=10000 <nb_ele>=8192 <nb_heterog>=0 <loading_angle>=90 ]" << std::endl;

  std::string sim_name = "Single material homogeneous interface controled rupture speed";

  // Geometry description
  UInt nb_time_steps = 10000; 
  UInt nb_elements = 8192;
  Real nu =  0.35;
  Real E = 5.3e9;
  Real cs = 1263;

  // Cut of the loaded material kernels
  UInt tcut = 100;
  
  // Loading case
  Real load = 3e6;
  Real psi = 90;
  Real phi = 0;

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

  UInt nb_heterog = 0;

  if(argc > 7)
    nb_heterog = std::atoi(argv[7]);

  if(argc > 8)
    psi = std::atof(argv[8]);

  std::cout << "./cst_speed_fracture " 
	    << "v_imposed:" << cr_speed << " "
	    << "load_file:" << load_file << " "
	    << "output folder:" << output_folder << " " 
	    << "nb_time steps:" << nb_time_steps << " " 
	    << "nb_elements:" << nb_elements << " "
	    << "nb_heterog:" << nb_heterog << " "
	    << "loading angle:" << psi << " "
	    << std::endl;
  
  // Cohesiv paramters
  Real crit_n_open = 0.02e-3;
  Real crit_s_open = 0.02e-3;
  Real max_s_str = 9e6;
  Real max_n_str = 9e6;
  Real wk_max_s_str = 4e6;
  Real wk_max_n_str = 4e6;
  Real sg_max_n_str = 14e6;
  Real sg_max_s_str = 14e6;

  Real G_length = crit_s_open*max_s_str/(load*load*M_PI)*E/(1-nu*nu);

  std::cout << "Griffith length: " << G_length << std::endl; 

  Real dom_size = 50*G_length;
  Real dx = dom_size/double(nb_elements);
  Real wall_position = 0.3*dom_size;
  UInt propagation_domain = 0.25*nb_elements;
  
  // Friction paramters
  Real regularized_time_scale = 0.1;
  Real coef_frict = 0.25;

  std::shared_ptr<ContactLaw> contactlaw = std::make_shared<RegularizedCoulombLaw>(coef_frict, regularized_time_scale, nb_elements);

  /* -------------------------------------------------------------------------- */

  SpectralModel model({nb_elements,1}, nb_time_steps, {dom_size,0.}, 
		      nu, nu, E, E, cs, cs, tcut, tcut, sim_name, output_folder); 

  // SimulationDriver object helping to launch controlled-speed simulation
  SimulationDriver sim_driver(model, cr_speed, 0.0);
  
  Interfacer<_linear_coupled_cohesive> interfacer(model);

  DataRegister::registerParameter("critical_normal_opening",crit_n_open);
  DataRegister::registerParameter("critical_shear_opening",crit_s_open);
  DataRegister::registerParameter("max_normal_strength",max_n_str);
  DataRegister::registerParameter("max_shear_strength",max_s_str);
  interfacer.createUniformInterface();
  interfacer.createThroughCrack(0.,5*dx);
  
  if(nb_heterog>0) {
    
    interfacer.createThroughMultiPolAsperity(5*G_length, 11.25*G_length, nb_heterog,
					     (sg_max_s_str-max_s_str)/max_s_str, 
					     (sg_max_n_str-max_n_str)/max_n_str,
					     0.,0.,true);
  }
  
  interfacer.createThroughWall(wall_position,dom_size);

  std::shared_ptr<CohesiveLaw> cohesive_law = std::dynamic_pointer_cast<CohesiveLaw>(model.getInterfaceLaw());
  cohesive_law->preventSurfaceOverlapping(contactlaw);

  if(write){
    // Init algorithm to tailor constant speed loading conditions
    //sim_driver.initConstantSpeed(load, psi, phi, max_s_str);
    sim_driver.initConstantSpeed(load, psi, phi, max_s_str, 0.0, _space_control);
  }
  else {
    // Read an existing file to set loading conditions
    sim_driver.initLoadingFromFile(output_folder+load_file, _space_control, load, psi, phi);
  }
    
  DataDumper dumper(model);

  UInt integ_width = G_length/dx;
  
  dumper.initIntegratorsDumper("Energies.cra", {5*G_length,1.},{10*G_length,1.});
  dumper.initSurfingIntegratorsDumper("surfing_eq.cra", integ_width, 0, nb_elements,
				      {_radiated_energy},{"surfing_Eq"});
    
  if (write){
    dumper.initDumper("ST_Diagram_id.cra", _id_crack, 0.5, 1, 0, _binary);
    dumper.initDumper("ST_Diagram_shear_velo_jump.cra", _shear_velocity_jumps, 0.125, 4, 0, _binary);
  } else {
    dumper.initDumper("ST_Diagram_id.cra", _id_crack, 0.125, 4, 0, _binary);
    dumper.initDumper("ST_Diagram_shear_velo_jump.cra", _shear_velocity_jumps, 0.5, 1, 0, _binary);
  }

  UInt x_tip=0;

  sim_driver.launchCrack(0.,G_length,0.15);

  UInt t = 0;

  UInt max_t_step = 1.5*propagation_domain/(model.getBeta()*cr_speed);

  UInt print_info = 0.05*max_t_step;

  while ((x_tip<(propagation_domain))&&(t<max_t_step)) {

    // High level method embedding the resolution of one time step by the SpectralModel
    x_tip = sim_driver.solveStep();  

    dumper.dumpAll();

    if(t%100==0){
      //dumper.dump("ST_Diagram_id.cra");
    }

    if (t%((UInt)print_info)==0){
      std::cout << "Crack tip at " << (Real)(x_tip)/(Real)nb_elements*100 << "% "
		<< " or double max time " << (Real)t/(Real)max_t_step*100 << "% "
		<< std::endl;
    }
    ++t;
  }

  // Do not forget to create the resulting loading file
  if (write)
    sim_driver.writeLoading(output_folder+load_file);
  
  return 0;
}
