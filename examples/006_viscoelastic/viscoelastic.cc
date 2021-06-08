/**
 * @file   viscoelastic.cc
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Thu Nov 19 9:57:04 2020
 *
 * @brief  Viscoelastic cohesive law in 2d. 
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
#include "cohesive_law_viscoelastic.hh"
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

  std::cout << "./viscoelastic [ <output_folder_name>='./']" 
	    << std::endl;
  
  std::string sim_name = "Viscoelastic 2D";

  std::string output_folder = "output_viscoelastic/";
  
  std::string input_file = "input_viscoelastic.cra";

  DataRegister::readInputFile(input_file);

  // Geometry description 
  Real nu = DataRegister::getParameter<Real>("nu");
  Real E = DataRegister::getParameter<Real>("E");
  Real cs = DataRegister::getParameter<Real>("cs");

  // Cut of the loaded material kernels
  UInt tcut = DataRegister::getParameter<UInt>("tcut");
  
  // Loading case
  Real psi = 0;
  Real phi = 0;
  
  // Cohesive parameters
  Real crit_n_open = DataRegister::getParameter<Real>("critical_normal_opening");
  Real crit_s_open = DataRegister::getParameter<Real>("critical_shear_opening");
  Real max_s_str = DataRegister::getParameter<Real>("max_normal_strength");
  Real max_n_str = DataRegister::getParameter<Real>("max_shear_strength");
  Real lim_velocity = DataRegister::getParameter<Real>("lim_velocity");
  
  Real load = DataRegister::getParameter<Real>("load");
  
  //Real G_length = crit_n_open*max_n_str/(load*load*M_PI)*E/(1-nu*nu);

  Real dom_size = DataRegister::getParameter<Real>("dom_size");

  std::cout << "Summary of the current simulation: " << std::endl
	    << "Domain size: " << dom_size << std::endl; 

  UInt nb_time_steps = 0;
  UInt nb_elements = DataRegister::getParameter<UInt>("nb_elements");
  Real dx = dom_size/nb_elements;
  Real propagation_domain = 0.9*dom_size;
  
  // Friction paramters
  Real regularized_time_scale = 0.1;
  Real coef_frict = 0.25;
  std::shared_ptr<ContactLaw> contactlaw = std::make_shared<RegularizedCoulombLaw>(coef_frict, regularized_time_scale, nb_elements);
  
  std::cout << "Input parameters summary : "<< std::endl
	    << "Output directory : " << output_folder << std::endl
	    << "Loading : " << load << std::endl;

  /* -------------------------------------------------------------------------- */
  
  SpectralModel model(nb_elements, nb_time_steps, dom_size, nu,
		      E, cs, tcut, sim_name, output_folder);
  
  Real beta=0.1; //Stable time step coefficient
   
  model.initModel(beta);
  model.setLoadingCase(0, psi, phi);  

  Interfacer<_viscoelastic_coupled_cohesive> interfacer(model);
  //Interfacer<_linear_coupled_cohesive> interfacer(model);

  interfacer.createUniformInterface();

  CohesiveLawViscoelastic& cohesive_law = dynamic_cast<CohesiveLawViscoelastic&>((model.getInterfaceLaw()));
  //CohesiveLaw& cohesive_law = dynamic_cast<CohesiveLaw&>((model.getInterfaceLaw()));
  cohesive_law.preventSurfaceOverlapping(contactlaw);

  cohesive_law.initLinearFormulation();

  interfacer.createThroughCrack(0, 150*dx);
  interfacer.createThroughWall(propagation_domain, dom_size);
  
  model.updateLoads();
  model.initInterfaceFields();
  
  DataDumper dumper(model);

  dumper.initDumper("ST_Diagram_id.cra", _id_crack);
  dumper.initVectorDumper("ST_Diagram_top_y_displ.cra", _top_displacements, _y, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_x_displ.cra", _top_displacements, _x, 1.0, 1, 0, _text);
  dumper.initDumper("ST_Diagram_normal_velo_jump.cra", _normal_velocity_jumps);
  dumper.initDumper("ST_Diagram_normal_strength.cra", _normal_strength);
  dumper.initVectorDumper("ST_Diagram_normal_stress.cra", _interface_tractions, _y, 1.0, 1, 0, _text);

  UInt propag_int = (UInt)(0.95*propagation_domain*nb_elements/dom_size);

  UInt x_tip = 0;
  UInt t = 0;

  UInt t_max = 20000;
    
  while (t < t_max) {
    
    UInt top = 8000;
    
    UInt load_factor = std::min(t,top);
    
    model.setLoadingCase(load_factor*load, psi, phi);
    model.updateLoads();
    
    model.updateDisplacements();
    model.fftOnDisplacements();
    model.computeStress();
    model.computeInterfaceFields();
    model.increaseTimeStep();  
    
    if(t%100==0){
      std::cout << t << std::endl;
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
