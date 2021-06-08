/**
 * @file   mode_II_rate_and_state.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Sun Apr 29 10:40:06 2018
 *
 * @brief  
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
 */
/* -------------------------------------------------------------------------- */
#include "spectral_model.hh"
#include "simulation_driver.hh"
#include "interfacer.hh"
#include "cohesive_law.hh"
#include "rate_and_state_law.hh"
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

  std::cout << "./mode_II_rate_and_state [ <output_folder_name>='./' <load xtau_m> <dom_size>=0.5 <nb_elements>=512 <epsilon>=1e-7 ]"
	    << std::endl;

  std::string sim_name = "mode-II-rate-and-state";
  std::string output_folder = "./";
  Real load = 0.35e6;
  Real dom_size = 0.5;
  UInt nb_elements = 512;
  Real epsilon = 1e-7;

  if(argc > 1)
    output_folder = argv[1];
  //if(argc > 2)
    //load = std::atof(argv[2])*0.309665e6;  
    //load = std::atof(argv[2])*0.34475668e6;  
    //load = std::atof(argv[2])*1e3;  
  if(argc > 3)
    dom_size = std::atof(argv[3]);  
  if(argc > 4)
    nb_elements = std::atof(argv[4]);  
  if(argc > 5)
    epsilon = std::atof(argv[5]);  

  // Geometry description 
  Real mu = 3e9;
  Real rho = 1200;
  Real nu =  0.33;
  Real E = 2*mu*(1+nu);
  Real cs = sqrt(mu/rho);
  // Cut of the loaded material kernels
  UInt tcut = 100;
  
  // Loading case
  Real psi = 90;
  Real phi = 0;
      
  /* -------------------------------------------------------------------------- */

  SpectralModel model(nb_elements, 0, dom_size, nu, 
		      E, cs, tcut, sim_name, output_folder);
  
  model.initModel();
  model.readInputFile("input_PMMA.dat");

  Real theta = 5./18.;
  DataRegister::registerParameter("theta",theta);
  
  Real k = DataRegister::getParameter<Real>("k");
  //Real epsilon = DataRegister::getParameter("epsilon");
  Real v_min_fric = DataRegister::getParameter<Real>("min_ss_v");
  Real v_predictor = DataRegister::getParameter<Real>("v0_predictor");
  UInt t_char = DataRegister::getParameter<UInt>("dumping_every_t");

  std::cout << "Input parameters summary : "<< std::endl
	    << "Output directory : " << output_folder << std::endl
	    << "k : " <<  k << std::endl  
	    << "epsilon : " << epsilon << std::endl    
	    << "minimum-friction velocity : " << v_min_fric << std::endl
	    << "implicit velocity predictor : " << v_predictor << std::endl
	    << "dumping characteristic time : " << t_char << std::endl;

  DataRegister::out_parameters << "k " << k << std::endl;
  DataRegister::out_parameters << "epsilon " << epsilon << std::endl;

  Interfacer<_rate_and_state> interfacer(model);

  RateAndStateLaw& r_and_s = dynamic_cast<RateAndStateLaw&>((model.getInterfaceLaw()));
  
  r_and_s.initStateEvolution();
  interfacer.createUniformInterface();
  
  model.setLoadingCase(load, psi, phi);
  model.updateLoads();
    
  DataDumper dumper(model);

  dumper.initVectorDumper("ST_Diagram_top_z_velo.cra", _top_velocities, 0, 1.0, 1, 0, _binary);
  dumper.initVectorDumper("ST_Diagram_top_z_displ.cra", _top_displacements, 0, 1.0, 1, 0, _binary);
  dumper.initVectorDumper("ST_Diagram_shear_stress.cra", _interface_tractions, 0, 1.0, 1, 0, _binary);
  dumper.initDumper("ST_Diagram_state_variable.cra", _state_variable, 1.0, 1, 0, _binary);
  dumper.initDumper("ST_Diagram_fric_coef.cra", _friction_coefficient, 1.0, 1, 0, _binary);
  
  UInt t = 0;

  const CrackProfile * shear_velo_jump = model.readData(_shear_velocity_jumps);
  const std::vector<Real> * state = model.readData(_state_variable);

  r_and_s.setVelocityPredictor({v_predictor,0.,0.});
  
  Real v_max=0.;
  Real v_av,v_25,state_75,state_ss; 
  bool dynamic = false;

  while (v_25*cs<v_min_fric) {
    //while (t<50000) {

    model.updateDisplacements();
    model.fftOnDisplacements();
    model.computeStress();

    if (t==0) {
      model.initInterfaceFields();
      state_ss = (*state)[(int)(0.75*nb_elements)];
    }
    else
      model.computeInterfaceFields();

    state_75 = (*state)[(int)(0.75*nb_elements)];
    model.increaseTimeStep();

    /*if ((!dynamic)&&(state_75<0.2*state_ss)){
      dynamic = true;
      std::cout << "Dynamic rupture just started..." << std::endl;
      }*/

    if (t==5)
      r_and_s.perturbState(epsilon,k);

    if ((t%t_char==0)||((dynamic)&&(t%5==0)))
      dumper.dumpAll();
    
    if (t%t_char==0){

      v_av = (*shear_velo_jump)[0];
      v_25 = (*shear_velo_jump)[(int)(0.25*nb_elements)];
      v_max = shear_velo_jump->getMaxValue();

      std::cout << "Simulation at t " << t << " = " << model.getTime()*dom_size/cs << " [sec]"<< std::endl
		<< "sliding velocity at the observation point: " << v_25*cs << std::endl
		<< "state variable at the expected nucleation point: " << state_75/state_ss << std::endl
		<< "-> v_max-v_av = " << (v_max-v_av)*cs << std::endl;
      }
    ++t;
  }
   
  return 0;
}
