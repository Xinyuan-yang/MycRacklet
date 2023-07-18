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
  std::cout << "cs = " << cs << std::endl;
  // Cut of the loaded material kernels
  UInt tcut = 100; 
  
  // Loading case
  Real psi = 90.0;
  Real phi = 90.0;
  Real start = 0.0; // valeur de départ
  Real end = 2.25e6; // valeur finale
  int n = 5000; // nombre d'éléments souhaité
  // Créer un vecteur de n éléments initialisé à 0
  std::vector<double> loads(n, 0.0);
  // Calculer l'incrément entre les valeurs
  double increment = (end - start) / (n - 1);
  // Remplir le vecteur avec les valeurs désirées
  std::transform(loads.begin(), loads.end(), loads.begin(), [start, increment](double) mutable {
      double value = start;
      start += increment;
      return value;
  });

  // Cohesive parameters
  Real crit_n_open = 0.00025; //50.0e-5;
  Real crit_s_open = 0.00025; //50.0e-5;
  Real max_n_str = 2.5e6;
  Real max_s_str = 2.5e6;
  Real res_n_str = 0.2e6;
  Real res_s_str = 0.2e6;
  Real nor_op_factor = 0.2;
  Real shr_op_factor = 0.2;
  Real nor_str_factor = 6;
  Real shr_str_factor = 6;

  //Real Gc = 0.5*crit_s_open*shr_op_factor*(max_s_str-res_s_str*shr_str_factor) + 0.5*crit_s_open*(1-shr_op_factor)*res_s_str*(shr_str_factor-1) + crit_s_open*res_s_str*(shr_str_factor-1);
  //Gc = 237500;
  //Real Gc = shr_op_factor*crit_s_open*(0.5*(max_s_str - shr_str_factor*res_s_str) + shr_str_factor*res_s_str - res_s_str) + 0.5*crit_s_open*(1-shr_op_factor)*(shr_str_factor*res_s_str-res_s_str);
  
  std::vector<double> op_list = {0.1, 0.5};
  //std::vector<double> str_list = {3.4e6, 3.2e6, 1.8e6, 1.6e6};
  std::vector<double> str_list = {2.0e6, 0.6e6};
  str_list.insert(str_list.begin(), max_s_str);
  str_list.insert(str_list.end(), res_s_str);
  op_list.insert(op_list.begin(), 0.0);
  op_list.insert(op_list.end(), 1.0);
  Real Gc = 0;
  for (int i=1; i<op_list.size(); i++) {
    Gc = Gc + crit_s_open*(op_list[i] - op_list[i-1])*(0.5*(str_list[i-1]-str_list[i]) + (str_list[i] - str_list.back()));
  }
  
  
  
  std::cout << "Gc =" << Gc << std::endl;
  op_list = {0.1, 0.5};
  //str_list = {3.4e6, 3.2e6, 1.8e6, 1.6e6};
  str_list = {2.0e6, 0.6e6};


  //Real G_length = 2*mu*crit_n_open*(max_n_str-res_n_str)/((load-res_n_str)*(load-res_n_str)*M_PI);
  Real G_length = 4*mu*Gc/(M_PI*std::pow(loads.back()-res_s_str, 2));
  //Real G_length = 4*mu*Gc/(M_PI*std::pow(loads.back()-0.25e6, 2));

    std::cout << "G_length =" << G_length << std::endl;
  
  Real dom_sizex = 15*G_length;
  Real dx = dom_sizex/(Real)(nex);

  //Real crack_size = 2*dx;
  Real crack_size = 2*G_length;
   
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
  SimulationDriver sim_driver(*model);

  model->setLoadingCase(0, psi, phi);

  Interfacer<_coupled_cohesive> interfacer(*model);   

  DataRegister::registerParameter("critical_normal_opening",crit_n_open);
  DataRegister::registerParameter("critical_shear_opening",crit_s_open);
  DataRegister::registerParameter("max_normal_strength",max_n_str);
  DataRegister::registerParameter("max_shear_strength",max_s_str);
  DataRegister::registerParameter("res_shear_strength",res_s_str);
  DataRegister::registerParameter("res_normal_strength",res_n_str);
  interfacer.createUniformInterface();

  interfacer.createThroughCrack((dom_sizex-crack_size)/2.,(dom_sizex+crack_size)/2.);    

  CohesiveLawAll& cohesive_law = dynamic_cast<CohesiveLawAll&>((model->getInterfaceLaw()));
  
  cohesive_law.preventSurfaceOverlapping(NULL);

  //cohesive_law.initRegularFormulation();
  //cohesive_law.initDualFormulation(nor_op_factor, shr_op_factor, nor_str_factor, shr_str_factor);  
  //cohesive_law.initTanhFormulation(0.5,0.15);
  cohesive_law.initMultiFormulation(op_list, str_list);

  //sim_driver.initConstantLoading(load, psi, phi);
    
  /* -------------------------------------------------------------------------- */
  //Set-up simulation  outputs
     
  //DataDumper dumper(*model);
  DataDumper dumper(*model);

  dumper.initVectorDumper("ST_Diagram_top_z_velo.cra", _top_velocities, 2, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_z_displ.cra", _top_displacements, 2, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_x_velo.cra", _top_velocities, 1, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_x_displ.cra", _top_displacements, 1, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_shear_stress.cra", _interface_tractions, 2, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_loading.cra", _top_loading, 2, 1.0, 1, 0, _text);
  dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1.0, 1, 0, _text);
  dumper.initDumper("ST_Diagram_normal_stress.cra", _interface_tractions, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_maximum_shear_strength.cra", _maximum_shear_strength, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_maximum_normal_strength.cra", _maximum_normal_strength, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_shear_vel.cra", _shear_velocity_jumps, 1.0, 1, 0);
  dumper.initVectorDumper("ST_Diagram_shear_tra.cra", _interface_tractions, 2, 1.0, 1, 0);
  dumper.initVectorDumper("ST_Diagram_normal_tra.cra", _interface_tractions, 1, 1.0, 1, 0);

  /* -------------------------------------------------------------------------- */
    
  //sim_driver.launchCrack(dom_sizex/2.,45*G_length,0.075,false);
  //sim_driver.launchCrack(dom_sizex/2.,1.75*G_length,0.075,false);

  model->updateLoads();
  model->initInterfaceFields();
  UInt i = 0;

  //while ((t < nb_time_steps)&&(x_tip<0.9*nex)) {
  while ((t < n)&&(x_tip<0.9*nex)) {

    //sim_driver.solveStep();
    model->setLoadingCase(loads[i], psi, phi);
    model->updateLoads();
    
    model->updateDisplacements();
    model->fftOnDisplacements();
    model->computeStress();
    model->computeInterfaceFields();
    model->increaseTimeStep();  
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
  ++i;
  }

  model->setLoadingCase(loads.back(), psi, phi);
  model->updateLoads();
  while ((t < nb_time_steps)&&(x_tip<0.9*nex)) {   
    model->updateDisplacements();
    model->fftOnDisplacements();
    model->computeStress();
    model->computeInterfaceFields();
    model->increaseTimeStep();  
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
  ++i;
  }
  //delete model;
  return 0;
}
