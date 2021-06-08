/**
 * @file   andrews_transition.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Aug 11 10:19:07 2016
 *
 * @brief  Study super-shear transition following Andrews(1976) formalism.
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

  std::cout << "./andrews_transition [ <output_folder_name>='./' <is_heterogeneous>=true <loading(x1e5)>=20.0 <nb_heterog=300> <toughness_ratio=5.0> ]" 
	    << std::endl;
  
  std::string sim_name = "Super-shear transition study";

  std::string output_folder = "./";
  bool heterog = true;
  Real load = 20.0;

  if(argc > 1)
    output_folder = argv[1];
  if(argc > 2) 
    if(std::atoi(argv[2])==0)
      heterog = false;
  if(argc > 3)
    load = std::atof(argv[3]);

  load*=1e5;

  // Geometry description 
  Real nu =  0.35;
  Real E = 5.3e9;
  Real cs = 1263;

  // Cut of the loaded material kernels
  UInt tcut = 100;
  
  // Loading case
  Real psi = 90;
  Real phi = 0;
   
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

  Real G_length = wk_crit_s_open*mean_max_s_str/(load*load*M_PI)*E/(1-nu*nu);
  Real S_ratio = (mean_max_s_str-load)/load;
  Real dom_size = 100*G_length;
  std::cout << "Summary of the current simulation: " << std::endl
	    << "Seismic ratio: " << S_ratio << std::endl 
	    << "Griffith length: " << G_length << std::endl; 

  //int nb_time_steps = 160000/4*dom_size;
  UInt nb_time_steps = 0;
  UInt nb_elements = 16384;
  Real dx = dom_size/nb_elements;
  Real propagation_domain = 0.9*dom_size;
  Real crack_size = 0.6*G_length;
 
  // Number of heterogeneous bands
  UInt nb_heterog = 300;
  Real toughness_ratio = 5.0;

  if(argc > 4)
    nb_heterog = std::atoi(argv[4]);
  if(argc > 5)
    toughness_ratio = std::atof(argv[5]);

  sg_max_s_str = 2*mean_max_s_str*toughness_ratio/(1+toughness_ratio);
  sg_max_n_str = 2*mean_max_n_str*toughness_ratio/(1+toughness_ratio);

  Real crk_srt = crack_size;
 
  // Friction paramters
  Real regularized_time_scale = 0.1;
  Real coef_frict = 0.25;
  std::shared_ptr<ContactLaw> contactlaw = std::make_shared<RegularizedCoulombLaw>(coef_frict, regularized_time_scale, nb_elements);
   
  // Output point position
  UInt nb_obs_points = 25;
  UInt start = (UInt)(crk_srt/dom_size*nb_elements);
  UInt end = (UInt)(propagation_domain/dom_size*nb_elements);

  std::cout << "Input parameters summary : "<< std::endl
	    << "Output directory : " << output_folder << std::endl
	    << "Is heterogeneous : " << heterog << std::endl
	    << "Loading : " << load << std::endl
	    << "Number of heterogeneities : " << nb_heterog << std::endl
	    << "Toughness ratio : " << toughness_ratio << std::endl;

  /* -------------------------------------------------------------------------- */
  
  SpectralModel model(nb_elements, nb_time_steps, dom_size, nu,
		      E, cs, tcut, sim_name, output_folder);
  
  SimulationDriver sim_driver(model);

  Interfacer<_linear_coupled_cohesive> interfacer(model);

  DataRegister::registerParameter("critical_normal_opening",wk_crit_n_open);
  DataRegister::registerParameter("critical_shear_opening",wk_crit_s_open);
  DataRegister::registerParameter("max_normal_strength",mean_max_n_str);
  DataRegister::registerParameter("max_shear_strength",mean_max_s_str);
  interfacer.createUniformInterface();

  interfacer.createThroughCrack(0., 5*dx);
  UInt x_end;

  if(heterog) {
  x_end = interfacer.createThroughMultiPolAsperity(2*G_length, propagation_domain, nb_heterog,
							(sg_max_s_str-mean_max_s_str)/mean_max_n_str, 
							(sg_max_s_str-mean_max_n_str)/mean_max_s_str,
							0.,0.,true);
  interfacer.createThroughWall(x_end*dx,dom_size);  
  }
  else 
    interfacer.createThroughWall(propagation_domain,dom_size);  

  CohesiveLaw& cohesive_law = dynamic_cast<CohesiveLaw&>((model.getInterfaceLaw()));
  cohesive_law.preventSurfaceOverlapping(contactlaw);

  sim_driver.initConstantLoading(load,psi,phi);
  
  DataDumper dumper(model);

  dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1.0, 1, 0, _binary);
  dumper.initDumper("ST_Diagram_shear_velo_jump.cra", _shear_velocity_jumps, 1.0, 1, 0, _binary);
  dumper.initDumper("ST_Diagram_shear_strength.cra", _shear_strength, 1.0, 1, 0, _binary);
  dumper.initPointsDumper("Point_history.dat", start, end, nb_obs_points); 

  UInt propag_int = (UInt)(0.95*propagation_domain*nb_elements/dom_size);

  UInt x_tip = 0;
  UInt t = 0;

  sim_driver.launchCrack(0.,G_length,0.2);
  
  while (x_tip < propag_int) {

    sim_driver.solveStep();
  
    if(t%5==0){
      dumper.dumpAll();
    }
    
    x_tip = model.getCrackTipPosition(0.,nb_elements);
    if (x_tip%((UInt)(0.05*nb_elements))==0){
      std::cout << "Crack position at " << x_tip*(Real)(dx) << " over " << propagation_domain<< std::endl;
    }
    ++t;
  }
   
  return 0;
}
