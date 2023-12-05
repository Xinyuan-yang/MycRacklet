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
 *
 */

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
#include <iomanip>
#include <sys/stat.h>
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]){

  // Note : Construct the pre-integrated material kernels before running this simulation
  // Use "invert_serial.f" to construct kernel files
  
  std::string output_folder=argv[1];
  // Geometry description

  UInt nb_time_steps = std::atoi(argv[3]); 
  UInt nex = std::atoi(argv[2]); 
  Real mu = 3e9;
  Real rho = 1200;
  Real nu =  0.33;
  Real E = 2*mu*(1+nu);
  Real cs = sqrt(mu/rho);
  std::cout << "cs = " << cs << std::endl;
  // Cut of the loaded material kernels
  UInt tcut = 100; 
  
  // Loading case
  Real load = 2.5e6; //2.25e6;
  Real psi = 90.0;
  Real phi = 90.0;

  Real sigma_0 = 1e7;

  // Cohesive parameters
  Real coeff_s = 0.25;
  Real coeff_d = 0.025;
  Real coeff_int = 0.2;
  Real crit_n_open = 0.0005;
  Real crit_s_open = 0.0005;
  Real crit_int_open = 0.00001;
  //Real nor_op_factor = 0.02;
  //Real shr_op_factor = 0.02;
  //Real nor_str_factor = 8;
  //Real shr_str_factor = 8;

  Real Gc = 0.5*crit_s_open*(coeff_s - coeff_d)*sigma_0;

  //Real Gc = sigma_0 * (0.5* (coeff_s - coeff_int)*crit_int_open + (coeff_int-coeff_d)*crit_int_open + 0.5*(coeff_int-coeff_d)*(crit_s_open-crit_int_open));
  
  std::cout << "Gc =" << Gc << std::endl;

  Real G_length = 4*mu*Gc/(M_PI*std::pow(load-coeff_d*sigma_0, 2));

  std::cout << "G_length =" << G_length << std::endl;
  
  Real dom_sizex = 15*G_length;
  Real dx = dom_sizex/(Real)(nex);

  //Real crack_size = 2*dx;
  Real crack_size = 1*G_length;
   
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
  
  model = new SpectralModel(nex, 0, dom_sizex,
				nu, E, cs, tcut,
				sim_name, output_folder);      

  //Real beta=0.002;
  //SimulationDriver sim_driver(*model, beta=beta);
  //SimulationDriver sim_driver(*model);
  model->initModel();

  DataRegister::registerParameter("critical_normal_opening",crit_n_open);
  DataRegister::registerParameter("critical_shear_opening",crit_s_open);
  DataRegister::registerParameter("static_friction_coefficient",coeff_s);
  DataRegister::registerParameter("dynamic_friction_coefficient",coeff_d);
  DataRegister::registerParameter("uniform_contact_pressure",sigma_0);

  // Dual
  DataRegister::registerParameter("int_friction_coefficient",coeff_int);
  DataRegister::registerParameter("critical_int_opening",crit_int_open);

    /* -------------------------------------------------------------------------- */
  //Set-up simulation  outputs
     
  //DataDumper dumper(*model);
  DataDumper dumper(*model);
  dumper.initVectorDumper("ST_Diagram_top_z_velo.cra", _top_velocities, 2, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_z_displ.cra", _top_displacements, 2, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_x_velo.cra", _top_velocities, 1, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_x_displ.cra", _top_displacements, 1, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_shear_stress.cra", _interface_tractions, 2, 1.0, 1, 0, _text);
  dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1.0, 1, 0, _text);
  dumper.initDumper("ST_Diagram_normal_stress.cra", _interface_tractions, 1.0, 1, 0);
  //dumper.initDumper("ST_Diagram_maximum_shear_strength.cra", _maximum_shear_strength, 1.0, 1, 0);
  //dumper.initDumper("ST_Diagram_maximum_normal_strength.cra", _maximum_normal_strength, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_shear_vel.cra", _shear_velocity_jumps, 1.0, 1, 0);
  dumper.initVectorDumper("ST_Diagram_shear_tra.cra", _interface_tractions, 2, 1.0, 1, 0);
  dumper.initVectorDumper("ST_Diagram_normal_tra.cra", _interface_tractions, 1, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_shear_displ.cra", _shear_displacement_jumps, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_static_coefficient.cra", _static_friction_coefficient, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_dynamic_coefficient.cra", _dynamic_friction_coefficient, 1, 1.0, 1, 0, _text);
  dumper.initDumper("ST_Diagram_fric_coef.cra", _friction_coefficient, 1.0, 1, 0, _text);
  dumper.initDumper("ST_Diagram_frictional_strength.cra", _frictional_strength, 1.0, 1, 0, _text);

  //Interfacer<_cohesive_coulomb> interfacer(*model);   
  Interfacer<_dual_cohesive_coulomb> interfacer(*model);   
  CohesiveLawCoulomb& cohesive_law = dynamic_cast<CohesiveLawCoulomb&>((model->getInterfaceLaw()));
  interfacer.createUniformInterface();
  
  interfacer.createThroughCrack((dom_sizex-crack_size)/2.,(dom_sizex+crack_size)/2.);    
  //CohesiveLawAll& cohesive_law = dynamic_cast<CohesiveLawAll&>((model->getInterfaceLaw()));
  dumper.dumpAll();
  //cohesive_law.preventSurfaceOverlapping(NULL);

  //cohesive_law.initRegularFormulation();
  //cohesive_law.initDualFormulation(nor_op_factor, shr_op_factor, nor_str_factor, shr_str_factor);  
  //cohesive_law.initTanhFormulation(0.5,0.15);
  //cohesive_law.initMultiFormulation(op_list, str_list);
  //sim_driver.initConstantLoading(load, psi, phi);
  model->setLoadingCase(load, psi, phi);
  model->updateLoads();
  model->initInterfaceFields();
  /* -------------------------------------------------------------------------- */
  //sim_driver.launchCrack(dom_sizex/2.,45*G_length,0.075,false);
  //sim_driver.launchCrack(dom_sizex/2.,1.75*G_length,0.075,false);
  while ((t < nb_time_steps)&&(x_tip<0.9*nex)) {
    //sim_driver.solveStep();
    model->updateDisplacements(); 
    model->fftOnDisplacements();
    model->computeStress();
    model->computeInterfaceFields();
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

    model->increaseTimeStep();

    ++t;
  }
  //delete model;
  return 0;
}
