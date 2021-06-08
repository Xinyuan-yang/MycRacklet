/**
 * @file   test_alu_homa3d_restart.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Oct 19 10:59:27 2017
 *
 * @brief  Testing restart between 2d and 3d simulations
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

int main(){

  // Note : Construct the pre-integrated material kernels before running this simulation
  // Use "invert_serial.f" to construct kernel files

  // Geometry description
  UInt nb_time_steps = 200; 
  UInt nex = 16;
  UInt nez = 8;
  Real dom_sizex = 1.0;
  Real dom_sizez = 1.0;
  Real crack_size=2*dom_sizex/nex;
  Real nu_mtl =  0.33;
  Real nu_poly = 0.35;
  Real E_mtl = 71e9;
  Real E_poly = 5.3e9;
  Real cs_mtl = 3100;
  Real cs_poly = 1263; 

  // Cut of the loaded material kernels
  UInt tcut_mtl = 100; //mis dans les kernels setup
  UInt tcut_poly = 100;
  
  // Loading case
  Real load = 3e6;//3e6
  Real psi = 75;
  Real phi = 45;
   
  // Cohesive paramters
  Real crit_n_open = 0.02e-3;
  Real crit_s_open = 0.02e-3;
  Real max_s_str = 5e6;
  Real max_n_str = 5e6;
   
  // Friction paramters
  Real regularized_time_scale = 0.1;
  Real coef_frict = 0.25;
  
  /* -------------------------------------------------------------------------- */

  mkdir("2d_outputs",0777);
  mkdir("3d_outputs",0777);
  
  UInt print = 0.01*nb_time_steps;
  
  for (UInt step = 0; step < 2; ++step) {
    
    std::shared_ptr<ContactLaw> contactlaw = std::make_shared<RegularizedCoulombLaw>(coef_frict, regularized_time_scale, nex*nez);
    SpectralModel * model;

    UInt t_start, t_end;
      
  
    if (step==0) {
      nez = 1;
      model = new SpectralModel(nex, nb_time_steps, dom_sizex,
				nu_mtl,nu_poly, E_mtl, E_poly, cs_mtl, cs_poly, 
				tcut_mtl, tcut_poly, 
				"test_alu_homa2d","2d_outputs/");
      t_start = 0;
      t_end = 101;
    }
    else {
      nez = 8;
      model = new SpectralModel({nex,nez}, nb_time_steps, {dom_sizex,dom_sizez},
				nu_mtl,nu_poly, E_mtl, E_poly, cs_mtl, cs_poly, 
				tcut_mtl, tcut_poly,
				"test_alu_homa3d","3d_outputs/");
      t_start = 101;
      t_end = nb_time_steps;
    }
    
    SimulationDriver driver(*model,0.37);
  
    Interfacer<_linear_coupled_cohesive> interfacer(*model);
    interfacer.createThroughCenteredCrack(crack_size, crit_n_open, crit_s_open, max_n_str, max_s_str);

    CohesiveLaw& cohesive_law = dynamic_cast<CohesiveLaw&>((model->getInterfaceLaw()));
    cohesive_law.preventSurfaceOverlapping(contactlaw);
    
    driver.initConstantLoading(load, psi, phi);
  
    const CrackProfile * t_displacements = model->readData(_top_displacements);
    const CrackProfile * b_displacements = model->readData(_bottom_displacements);
    const std::vector<Real> * nor_strength = model->readData(_normal_strength);
    const std::vector<Real> * shr_strength = model->readData(_shear_strength);

    std::vector<UInt> integ_points_left(nez*nex/2);
  
    UInt local_index = 0;
    for (UInt ix = 0; ix < nex/2; ++ix) {
      for (UInt iz = 0; iz < nez; ++iz) {
	UInt l_index = ix+iz*nex;
	integ_points_left[local_index] = l_index;
	++local_index;
      }
    }

    std::vector<Real> dx = model->getElementSize();
  
    Integrator E_s(integ_points_left, _shear_fracture_energy, 0., dx[0]*dx[1]);
    model->registerComputer("efrac_shear_left",&E_s);
    Integrator E_n(integ_points_left, _normal_fracture_energy, 0., dx[0]*dx[1]);
    model->registerComputer("efrac_normal_left",&E_n);
    Integrator E_fr(integ_points_left, _frictional_energy, 0., dx[0]*dx[1]);
    model->registerComputer("efric_left",&E_fr);

    UInt chkptx = 0.75*nex;
    UInt chkptz = 0;
    UInt wdth = 15;
    std::cout  << std::setw(wdth) << "delta_s" 
	       << std::setw(wdth) << "delta_n" 
	       << std::setw(wdth) << "tau_str^s" 
	       << std::setw(wdth) << "tau_str^n" 
	       << std::setw(wdth) << "E_s" 
	       << std::setw(wdth) << "E_n" 
	       << std::setw(wdth) << "E_fr" 
	       << std::endl;

    std::cout.precision(6);

    if(step==1) {
      DataRegister::restart_dir = "../2d_outputs/restart_files/";
      model->restartModel(true);
    }

    for (UInt t = t_start; t < t_end ; ++t) {

      driver.solveStep();
      
      if (print == (UInt)(0.1*nb_time_steps)) {
	std::cout << std::setw(wdth) << (*t_displacements)[(chkptx+nex*chkptz)*3]-(*b_displacements)[(chkptx+nex*chkptz)*3]  
		  << std::setw(wdth) << (*t_displacements)[(chkptx+nex*chkptz)*3+1]-(*b_displacements)[(chkptx+nex*chkptz)*3+1] 
		  << std::setw(wdth) << (*shr_strength)[(chkptx+nex*chkptz)] 
		  << std::setw(wdth) << (*nor_strength)[(chkptx+nex*chkptz)] 
		  << std::setw(wdth) << E_s.getIntegration() 
		  << std::setw(wdth) << E_n.getIntegration() 
		  << std::setw(wdth) << E_fr.getIntegration()
		  << std::endl;
       
	print=0;
      }
      ++print;
    }
    if(step==0)
      model->pauseModel();
    
    delete model;
  }
  return 0;
}

