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

/* -------------------------------------------------------------------------- */
#include "spectral_model.hh"
#include "simulation_driver.hh"
#include "interfacer.hh"
#include "cohesive_law.hh"
#include "cohesive_law_all.hh"
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
  Real mu = 3e9;
  Real rho = 1200;
  Real nu =  0.33;
  Real E = 2*mu*(1+nu);
  Real cs = sqrt(mu/rho);

  // Cut of the loaded material kernels
  UInt tcut = 100;
  
  // Loading case
  Real load = 4.5e6;
  Real psi = 90;
  Real phi = 90;

  // Target crack speed v=cr_speed*c_s
  Real cr_speed = std::atof(argv[1]);
  // True=simulation to tailor the loading file
  
  std::string output_folder = argv[4];
    
  UInt nb_time_steps = std::atoi(argv[3]);

  UInt nb_elements = std::atoi(argv[2]);

  std::cout << "./cst_speed_fracture " 
	    << "v_imposed:" << cr_speed << " "
	    << "output folder:" << output_folder << " " 
	    << "nb_time steps:" << nb_time_steps << " " 
	    << "nb_elements:" << nb_elements << " "
	    << "loading angle:" << psi << " "
	    << std::endl;
  
  // Cohesiv paramters
  Real crit_n_open = 50.0e-5;
  Real crit_s_open = 50.0e-5;
  Real max_n_str = 5e6;
  Real max_s_str = 5e6;
  Real res_n_str = 0.25e6;
  Real res_s_str = 0.25e6;
  Real nor_op_factor = 0.2;
  Real shr_op_factor = 0.2;
  Real nor_str_factor = 5;
  Real shr_str_factor = 5;

  Real Gc = 0.5*crit_s_open*shr_op_factor*(max_s_str-res_s_str*shr_str_factor) + 0.5*crit_s_open*(1-shr_op_factor)*res_s_str*(shr_str_factor-1) + crit_s_open*res_s_str*(shr_str_factor-1);
  std::cout << "Gc =" << Gc << std::endl;
  Real G_length = 4*mu*Gc/(M_PI*std::pow(load-res_s_str, 2));
  std::cout << "Griffith length: " << G_length << std::endl; 

  Real dom_size = 15*G_length;
  Real dx = dom_size/double(nb_elements);
  Real crack_size = 2*dx;

  UInt t = 0;
  UInt x_tip=0;
  UInt x_lap = 0.05*nb_elements;

  SpectralModel * model;
  
  model = new SpectralModel(nb_elements, 0, dom_size,
				nu, E, cs, tcut,
				sim_name, output_folder);      

  SimulationDriver sim_driver(*model, cr_speed, dom_size/2);
  //SimulationDriver sim_driver(*model);

  // SimulationDriver object helping to launch controlled-speed simulation
  
  Interfacer<_coupled_cohesive> interfacer(*model);   

  DataRegister::registerParameter("critical_normal_opening",crit_n_open);
  DataRegister::registerParameter("critical_shear_opening",crit_s_open);
  DataRegister::registerParameter("max_normal_strength",max_n_str);
  DataRegister::registerParameter("max_shear_strength",max_s_str);
  DataRegister::registerParameter("res_shear_strength",res_s_str);
  DataRegister::registerParameter("res_normal_strength",res_n_str);
  interfacer.createUniformInterface();
  interfacer.createThroughCrack((dom_size-crack_size)/2.,(dom_size+crack_size)/2.);    
  


  CohesiveLawAll& cohesive_law = dynamic_cast<CohesiveLawAll&>((model->getInterfaceLaw()));  
  cohesive_law.preventSurfaceOverlapping(NULL);

  cohesive_law.initRegularFormulation();
  //cohesive_law.initDualFormulation(nor_op_factor, shr_op_factor, nor_str_factor, shr_str_factor);

  // Init algorithm to tailor constant speed loading conditions
  //sim_driver.initConstantSpeed(load, psi, phi, max_s_str);
  //sim_driver.initConstantLoading(load, psi, phi);
  sim_driver.initConstantSpeed(load, psi, phi, max_s_str, 0.0, _space_control);

  DataDumper dumper(*model);

  dumper.initVectorDumper("ST_Diagram_top_z_velo.cra", _top_velocities, 2, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_z_displ.cra", _top_displacements, 2, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_x_velo.cra", _top_velocities, 1, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_x_displ.cra", _top_displacements, 1, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_shear_stress.cra", _interface_tractions, 2, 1.0, 1, 0, _text);
  dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1.0, 1, 0, _text);
  dumper.initDumper("ST_Diagram_normal_stress.cra", _interface_tractions, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_maximum_shear_strength.cra", _maximum_shear_strength, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_maximum_normal_strength.cra", _maximum_normal_strength, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_shear_vel.cra", _shear_velocity_jumps, 1.0, 1, 0);
  dumper.initVectorDumper("ST_Diagram_shear_tra.cra", _interface_tractions, 2, 1.0, 1, 0);
  dumper.initVectorDumper("ST_Diagram_normal_tra.cra", _interface_tractions, 1, 1.0, 1, 0);

  sim_driver.launchCrack(dom_size/2.,4*G_length,0.075,false);

  while ((t < nb_time_steps)&&(x_tip<0.9*nb_elements)) {

    sim_driver.solveStep();
    x_tip = model->getCrackTipPosition(nb_elements/2,nb_elements);

    if (t%10==0){
      dumper.dumpAll();
    }

    if ((x_tip>x_lap)||(t%(UInt)(0.05*nb_time_steps)==0)) {
      std::cout << "Process at " << (Real)t/(Real)nb_time_steps*100 << "% " << std::endl;
      std::cout << "Crack at " << 100*x_tip/(Real)(nb_elements) << "% " << std::endl;
      std::cout << std::endl;
      
      if (x_tip>x_lap)
	x_lap += 0.05*nb_elements;
    }

    ++t;
    
  }
  //delete model;
  return 0;
}
