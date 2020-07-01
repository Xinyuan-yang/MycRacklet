/**
 * @file   mode_III_slip_weakening.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Mar 19 15:47:26 2018
 *
 * @brief  Mode-III cract tip equation of motion
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
#include "coulomb_law.hh"
#include "regularized_coulomb_law.hh"
#include "data_dumper.hh"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <sys/stat.h>
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]){

  // Note : Construct the pre-integrated material kernels before running this simulation
  // Use "invert_serial.f" to construct kernel files

  std::cout << "./mode_III_slip_weakening <output_folder_name> <nb_ele_x> <nb_time_steps>" << std::endl;

  std::string output_folder=argv[1];
   
  // Geometry description
  UInt nb_time_steps = std::atoi(argv[3]); 
  UInt nex = std::atoi(argv[2]); 
  Real mu = 3e9;
  Real rho = 1200;
  Real nu =  0.33;
  Real E = 2*mu*(1+nu);
  Real cs = sqrt(mu/rho);
  // Cut of the loaded material kernels
  UInt tcut = 100; 
  
  // Loading case
  Real load = 1e6;
  Real psi = 90.0;
  Real phi = 90.0;

  // Cohesive parameters
  Real crit_n_open = 0.02e-3;
  Real crit_s_open = 0.02e-3;
  Real max_n_str = 5e6;
  Real max_s_str = 5e6;

  Real G_length = 2*mu*crit_n_open*max_n_str/(load*load*M_PI);
  
  Real dom_sizex = 15*G_length;
  Real dx = dom_sizex/(Real)(nex);

  Real crack_size = 2*dx;
   
  std::string sim_name = "Mode-III crack tip equation of motion";

  std::cout << "./mode_III_rate_and_state " 
	    << "output folder: " << output_folder << " " 
	    << "nb_elements alog x: " << nex << " "
	    << "nb_time_steps: " << nb_time_steps << " "
	    << "griffith crack length: " << G_length << " "
	    << std::endl;
   
  /* -------------------------------------------------------------------------- */

  UInt t = 0;
  UInt x_tip=0;
  UInt x_lap = 0.05*nex;

  SpectralModel * model;
    
  model = new SpectralModel({nex,1}, 0, {dom_sizex,0.},
				nu, nu, E, E, cs, cs, tcut, tcut,
				sim_name, output_folder);      
  
  SimulationDriver sim_driver(*model);
   
  Interfacer<_linear_coupled_cohesive> interfacer(*model);   

  DataRegister::registerParameter("critical_normal_opening",crit_n_open);
  DataRegister::registerParameter("critical_shear_opening",crit_s_open);
  DataRegister::registerParameter("max_normal_strength",max_n_str);
  DataRegister::registerParameter("max_shear_strength",max_s_str);
  interfacer.createUniformInterface();
  interfacer.createThroughCrack((dom_sizex-crack_size)/2.,(dom_sizex+crack_size)/2.);
    

  CohesiveLaw& cohesive_law = dynamic_cast<CohesiveLaw&>((model->getInterfaceLaw()));
  cohesive_law.preventSurfaceOverlapping(NULL);
    
  sim_driver.initConstantLoading(load, psi, phi);
    
  /* -------------------------------------------------------------------------- */
  //Set-up simulation  outputs
     
  DataDumper dumper(*model);

  dumper.initVectorDumper("ST_Diagram_top_z_velo.cra", _top_velocities, 2, 1.0, 1, 0, _binary);
  dumper.initVectorDumper("ST_Diagram_top_z_displ.cra", _top_displacements, 2, 1.0, 1, 0, _binary);
  dumper.initVectorDumper("ST_Diagram_shear_stress.cra", _interface_tractions, 2, 1.0, 1, 0, _binary);
  dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1.0, 1, 0, _binary);

  /* -------------------------------------------------------------------------- */
    
  sim_driver.launchCrack(dom_sizex/2.,1.75*G_length,0.075,false);

  while ((t < nb_time_steps)&&(x_tip<0.9*nex)) {

    sim_driver.solveStep();
    x_tip = model->getCrackTipPosition(nex/2,nex);
          
    if (t%10==0){
      dumper.dumpAll();
    }

    if ((x_tip>x_lap)||(t%(UInt)(0.05*nb_time_steps)==0)) {
      std::cout << "Process at " << (Real)t/(Real)nb_time_steps*100 << "% " << std::endl;
      std::cout << "Crack at " << 100*x_tip/(Real)(nex) << "% " << std::endl;
      std::cout << std::endl;
      
      if (x_tip>x_lap)
	x_lap += 0.05*nex;
    }

    ++t;
  }
  delete model;
  return 0;
}
