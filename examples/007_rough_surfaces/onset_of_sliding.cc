/**
 * @file   onset_of_sliding.cc
 * @author Fabian Barras <fabian.barras@mn.uio.no>
 * @date   Tue May 25 17:11:03 2021
 *
 * @brief  Nucleation of frictional rupture along two rough surfaces
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
 *
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
inline Real average(const std::vector<Real> * field) {
  Real sum = 0;
  Real nb_ele = 0;
  for (auto it = field->begin(); it != field->end(); ++it, ++nb_ele)
    sum += *it;
  
  return sum/nb_ele;
}

/* -------------------------------------------------------------------------- */
inline void createSineLoading(std::vector<Real> & loading, UInt nb_elem, Real min) {

  int half = nb_elem/2;
  
  for (int x = 0; x < nb_elem; ++x) {
    for (int z = 0; z < nb_elem; ++z) {
      Real radius = sqrt(Real(x-half)*Real(x-half)+Real(z-half)*Real(z-half));
      loading[x+z*nb_elem] = std::max(0.,cos(radius/nb_elem*M_PI))*(1.0-min) + min;
    }
  }
}

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]){

  std::cout << "./onset_of_sliding <output_folder_name> <nb_t_step>" << std::endl;
  
  std::string output_folder=argv[1];

  UInt nb_time_steps = std::atoi(argv[2]);
  
  // Geometry description
  // Nb elements
  UInt nex = 512;
  UInt nez = 512;

  //Microcontacts profile computed from elasto-plastic contact simulation (Tamaas)
  std::string roughness_profile_filename = "microcontacts_profile.txt";
  
  // Solid properties
  Real mu = 3e9;
  Real rho = 1200;
  Real nu =  0.33;
  Real E = 2*mu*(1+nu);
  Real cs = sqrt(mu/rho);
  // Cut-off used in the loaded material kernels
  UInt tcut = 100; 
 
  // Loading conditions
  // Initial fluid pressure
  Real tau0 = 1e6;
  Real tauf = 3e6;
  // Load increment spread over 80% of the time steps
  Real dtau = (tauf-tau0)/(0.9*nb_time_steps);
  // Initial far-field solid stress state
  Real psi = 90.0;
  Real phi = 0.0;

  // Cohesive parameters of the interface
  Real crit_n_open = 0.02e-3;
  Real crit_s_open = 0.02e-3;
  Real max_n_str = 5e6;
  Real max_s_str = 5e6;
  // The steps hereafter can be equivalently done using an input file
  DataRegister::registerParameter("critical_normal_opening",crit_n_open);
  DataRegister::registerParameter("critical_shear_opening",crit_s_open);
  DataRegister::registerParameter("max_normal_strength",max_n_str);
  DataRegister::registerParameter("max_shear_strength",max_s_str);

  // Grffith length
  Real G_length = crit_n_open*max_n_str/(tau0*tau0*M_PI)*E/(1-nu*nu);
  
  // Length of the domain (periodic boundary conditions at the domain edges)
  Real dom_size = 3*G_length;
  Real dx = dom_size/(Real)(nex);

  std::string sim_name = "Onset of frictional slip at the microcontacts scale";
  
  /* -------------------------------------------------------------------------- */
  // SpectralModel computing the elastodynamics
  SpectralModel * model; 
  model = new SpectralModel({nex,nez}, 0, {dom_size,dom_size},
				nu, nu, E, E, cs, cs, tcut, tcut,
				sim_name, output_folder);      

  // SimulationDriver handling the computation
  SimulationDriver sim_driver(*model);

  // Interfacer initialize the interface model and properties
  Interfacer<_linear_coupled_cohesive> interfacer(*model);
  interfacer.createUniformInterface();
  interfacer.insertPatternfromFile(roughness_profile_filename);
  
  // Set the far-field initial stress conditions in the solid
  std::vector<Real> & loading_ratio = model->getLoadingRatio();
  createSineLoading(loading_ratio, nex, 0.1);
  sim_driver.initConstantLoading(tau0, psi, phi);

  /* -------------------------------------------------------------------------- */
  //Set-up simulation  outputs
  // DataDumper handles output production 
  DataDumper dumper(*model);

  dumper.initDumper("ST_Diagram_shear_velocity_jump.cra", _shear_velocity_jumps, 1.0, 1, 0, _binary);
  dumper.initDumper("ST_Diagram_shear_strength.cra", _shear_strength, 1.0, 1, 0, _binary);
  const std::vector<Real> * strength = model->readData(_shear_strength);
  const CrackProfile * far_field_loading = model->readData(_top_loading);
  
  Real av_strength;
  Real av_loading;
  for(UInt t = 0; t < nb_time_steps; ++t) {

    // Solve one time step of the model
    sim_driver.solveStep();
    
    if (t%(UInt)(0.01*nb_time_steps)==0) {
      av_strength = average(strength);
      av_loading = 3*average(&(far_field_loading->getValues()));
      std::cout << "Process at " << (Real)t/(Real)nb_time_steps*100 << "% " << std::endl;
      std::cout << "Average frictional strength: " << av_strength
		<< "      Average shear loading: " << av_loading
		<< std::endl;
      std::cout << std::endl;
    }

    //Wait a couple of time steps before increasing the load
    if(t>0.1*nb_time_steps) {
      if(t%50==0)
	dumper.dumpAll();
      model->incrementLoad(dtau,0);
    }
    ++t;
  }
  
  delete model;

  return 0;
}
