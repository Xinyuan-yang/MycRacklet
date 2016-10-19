/**
 * @file   single_mat_heterog.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Sep 18 14:19:07 2014
 *
 * @brief  Heterogeneous striped interface between Homlalite bulks 
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

  std::cout << "./striped_interface [ <output_folder_name>='./' <is_heterogeneous>=true ]" 
	    << std::endl;

  
  std::string sim_name = "Stripes of heterogeneities meeting spontaneaous Mode-II crack";

  // Geometry description
  UInt nb_time_steps = 8000; 
  UInt nb_elements = 4096;
  Real crack_size = 0.1;
  Real nu =  0.35;
  Real E = 5.3e9;
  Real cs = 1263;

  // Cut of the loaded material kernels
  UInt tcut = 100;
  
  // Loading case
  Real load = 2e6;
  Real psi = 90;
  Real phi = 0;
  UInt l_index = 1;
   
  // Cohesive paramters
  Real wk_crit_n_open = 0.02e-3;
  Real wk_crit_s_open = 0.02e-3;
  Real wk_max_s_str = 3e6;
  Real wk_max_n_str = 3e6;
  Real sg_crit_n_open = wk_crit_s_open;
  Real sg_crit_s_open = wk_crit_n_open;
  Real sg_max_n_str = 15e6;
  Real sg_max_s_str = 15e6;
  Real mean_max_n_str=0.5*(sg_max_n_str+wk_max_n_str);
  Real mean_max_s_str=0.5*(sg_max_s_str+wk_max_s_str);

  // Critical stable crack size according to LEFM (=Griffith crack length)
  Real G_length = wk_crit_s_open*mean_max_s_str/(load*load*M_PI)*E/(1-nu*nu);
  Real dom_size = 51*G_length; //in parallel with andrews_transition study (98+2)*G_length
  Real propagation_domain = 0.9*dom_size;
  
  // Number of heterogeneous bands
  UInt nb_heterog = 50;
  
  FractureLaw * fracturelaw; 

  std::string output_folder = "./";
  bool heterog = true;

  if(argc > 1)
    output_folder = argv[1];
  if(argc > 2) 
    if(std::atoi(argv[2])==0)
      heterog = false;

  // Friction paramters
  bool overlap = 1;
  Real regularized_time_scale = 0.1;
  Real coef_frict = 0.25;
  ContactLaw * contactlaw = new RegularizedCoulombLaw(coef_frict, regularized_time_scale, nb_elements);

  /* -------------------------------------------------------------------------- */

  SpectralModel model({nb_elements,1}, nb_time_steps, {dom_size,0.}, 
		      nu, nu, E, E, cs, cs, tcut, tcut, overlap, l_index, 
		      fracturelaw, contactlaw, sim_name, output_folder);
  
  SimulationDriver sim_driver(model);
  
  Interfacer<_linear_coupled_cohesive> interfacer(model);
  interfacer.createUniformInterface(wk_crit_n_open, wk_crit_s_open, 
				    mean_max_n_str, mean_max_s_str);
  interfacer.createThroughCrack(0., 0.02*G_length);
  
  if(heterog)
    interfacer.createThroughMultiPolAsperity(2*G_length, propagation_domain, nb_heterog,
					     (sg_max_s_str-mean_max_s_str)/mean_max_n_str, 
					     (sg_max_s_str-mean_max_n_str)/mean_max_s_str,
					     0.,0.,true);

  interfacer.createThroughWall(propagation_domain,dom_size);  
  interfacer.applyInterfaceCreation();
  
  DataDumper dumper(model);

  dumper.initEnergetics("Energy.cra");
  std::string crack_id = "ST_Diagram_id.cra";
  std::string shr_v_jp = "ST_Diagram_shear_velo_jump.cra";
  std::string nor_d_jp = "ST_Diagram_normal_displ_jump.cra";
  dumper.initDumper(crack_id, _id_crack, 0.75, 1, 0);
  dumper.initDumper(shr_v_jp, _shear_velocity_jumps, 0.75, 1, 0, _binary);
  dumper.initDumper(nor_d_jp, _shear_strength, 0.75, 1, 0, _binary);

  sim_driver.launchCrack(0.,G_length,0.2);
  
  for (UInt t = 0; t < nb_time_steps ; ++t) {

    sim_driver.solveStep();
    
    if(t%5==0){
      dumper.printEnergetics();
      dumper.dumpAll();
    }

    if (t%((UInt)(0.1*nb_time_steps))==0){
      std::cout << "Process at " << (Real)t/(Real)nb_time_steps*100 << "% " << std::endl;      
    }
  }

  delete contactlaw;
   
  return 0;
}

