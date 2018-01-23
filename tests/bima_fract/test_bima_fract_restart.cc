/**
 * @file   test_bima_fract_restart.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Oct 18 15:29:59 2017
 *
 * @brief  Test restart on bimaterial fracture problem
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
 */
/* -------------------------------------------------------------------------- */
#include "spectral_model.hh"
#include "cohesive_law.hh"
#include "coulomb_law.hh"
#include "regularized_coulomb_law.hh"
#include "interfacer.hh"
#include "data_dumper.hh"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <iomanip>
/* -------------------------------------------------------------------------- */

int main(){

  // Note : Construct the pre-integrated material kernels before running this simulation
  // Use "invert_serial.f" to construct kernel files

  // Geometry description
  UInt nb_time_steps = 3500;
  UInt nb_elements = 2048;
  Real dom_size = 1.0;
  Real crack_size = 0.05;
  Real nu_mtl =  0.33;
  Real nu_poly = 0.35;
  Real E_mtl = 71e9;
  Real E_poly = 5.3e9;
  Real cs_mtl = 3100;
  Real cs_poly = 1263; 

  // Cut of the loaded material kernels
  UInt tcut_mtl = 100;
  UInt tcut_poly = 100;
  
  // Loading case
  Real load = 3e6;
  Real psi = 75;
  Real phi = 0;
   
  // Cohesive paramters
  Real crit_n_open = 0.02e-3;
  Real crit_s_open = 0.02e-3;
  Real max_s_str = 5e6;
  Real max_n_str = 5e6;

  // Friction paramters
  bool overlap = 0;
  Real regularized_time_scale = 0.1;
  Real coef_frict = 0.25;

  /* -------------------------------------------------------------------------- */

  UInt print = 0.01*nb_time_steps;
  
  for (UInt step = 0; step < 2; ++step) {
    
    FractureLaw * fracturelaw = NULL;    
    ContactLaw * contactlaw = new RegularizedCoulombLaw(coef_frict, regularized_time_scale, nb_elements);
    SpectralModel model({nb_elements,1}, nb_time_steps, {dom_size,0.}, 
			nu_mtl, nu_poly,E_mtl, E_poly, cs_mtl, cs_poly, 
			tcut_mtl, tcut_poly, overlap, fracturelaw,
			contactlaw,"test_bima_fract");
 
    model.initModel(0.4);
    model.setLoadingCase(load, psi, phi);

    Interfacer<_linear_coupled_cohesive> interfacer(model);
    interfacer.createThroughCenteredCrack(crack_size, crit_n_open, crit_s_open, max_n_str, max_s_str);
    interfacer.applyInterfaceCreation();
   
    model.updateLoads();
    model.computeInitialVelocities();

    const CrackProfile * t_displacements = model.readData(_top_displacements);
    const CrackProfile * b_displacements = model.readData(_bottom_displacements);
    const std::vector<Real> * nor_strength = model.readData(_normal_strength);
    const std::vector<Real> * shr_strength = model.readData(_shear_strength);

    std::vector<UInt> integ_points_left(nb_elements/2);
    std::vector<UInt> integ_points(nb_elements);
  
    for (UInt i = 0; i < nb_elements/2; ++i) {
      integ_points_left[i] = i;
      integ_points[2*i] = 2*i;
      integ_points[2*i+1] = 2*i+1;
    }

    std::vector<Real> dx = model.getElementSize();
  
    Integrator E_s(integ_points_left, _shear_fracture_energy, 0., dx[0]*dx[1]);
    model.registerComputer("efrac_shear_left",&E_s);
    Integrator E_n(integ_points_left, _normal_fracture_energy, 0., dx[0]*dx[1]);
    model.registerComputer("efrac_normal_left",&E_n);
    Integrator E_fr(integ_points_left, _frictional_energy, 0., dx[0]*dx[1]);
    model.registerComputer("efric_left",&E_fr);

    Integrator E_q(integ_points, _radiated_energy, 0., dx[0]*dx[1]);
    model.registerComputer("radiated_energy",&E_q);

    UInt chkpt = 0.75*nb_elements;
    UInt wdth = 15;
    std::cout  << std::setw(wdth) << "delta_s" 
	       << std::setw(wdth) << "delta_n" 
	       << std::setw(wdth) << "tau_str^s" 
	       << std::setw(wdth) << "tau_str^n" 
	       << std::setw(wdth) << "E_s" 
	       << std::setw(wdth) << "E_n" 
	       << std::setw(wdth) << "E_fr"
	       << std::setw(wdth) << "E_q" 
	       << std::endl;

    std::cout.precision(6);

    DataDumper dumper(model);
 
    dumper.initIntegratorsDumper("energetics.cra",integ_points_left,
				 {_shear_fracture_energy,_normal_fracture_energy,_frictional_energy},
				 {"Efrac shear left","Efrac normal left","Efric left"});

    UInt t_start, t_end;
    
    if(step==1) {
      model.restartModel();
      t_start = 2001;
      t_end = nb_time_steps;
    }
    else {
      t_start = 0;
      t_end = 2001;
    }
  
    for (UInt t = t_start; t < t_end ; ++t) {

      model.updateDisplacements();
      model.updateMaterialProp();
      model.fftOnDisplacements();
      model.computeStress();
      model.computeVelocities();
      model.increaseTimeStep();

      if (print == (UInt)(0.01*nb_time_steps)) {

	dumper.dumpAll();
	std::cout << std::setw(wdth) << (*t_displacements)[chkpt*3]-(*b_displacements)[chkpt*3]  
		  << std::setw(wdth) << (*t_displacements)[chkpt*3+1]-(*b_displacements)[chkpt*3+1] 
		  << std::setw(wdth) << (*shr_strength)[chkpt] 
		  << std::setw(wdth) << (*nor_strength)[chkpt] 
		  << std::setw(wdth) << E_s.getIntegration() 
		  << std::setw(wdth) << E_n.getIntegration() 
		  << std::setw(wdth) << E_fr.getIntegration()
		  << std::setw(wdth) << E_q.getIntegration()  
		  << std::endl;

	print=0;
      }

      ++print;
    }

    if(step==0)
      model.pauseModel();
  }   

  return 0;
}
