/**
 * @file   test_bima_fract.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Mar 11 11:19:07 2013
 *
 * @brief  Test of 2D bimaterial debonding
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
#include "cohesive_law.hh"
#include "coulomb_law.hh"
#include "regularized_coulomb_law.hh"
#include "interfacer.hh"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <iomanip>
/* -------------------------------------------------------------------------- */

 int main(){

   // Note : Construct the pre-integrated material kernels before running this simulation
   // Use "invert_serial.f" to construct kernel files

   // Geometry description
   int nb_time_steps = 3500;
   int nb_elements = 2048;
   double dom_size = 1.0;
   double crack_size = 0.05;
   double nu_mtl =  0.33;
   double nu_poly = 0.35;
   double E_mtl = 71e9;
   double E_poly = 5.3e9;
   double cs_mtl = 3100;
   double cs_poly = 1263; 

   // Cut of the loaded material kernels
   int tcut_mtl = 100;
   int tcut_poly = 100;
  
   // Loading case
   double load = 3e6;
   double psi = 75;
   double phi = 0;
   unsigned int l_index = 1;
   
   // Cohesive paramters
   double crit_n_open = 0.02e-3;
   double crit_s_open = 0.02e-3;
   double max_s_str = 5e6;
   double max_n_str = 5e6;
   FractureLaw * fracturelaw = new CohesiveLaw(crit_n_open, crit_s_open, max_n_str, max_s_str);

   // Friction paramters
   bool overlap = 0;
   double regularized_time_scale = 0.1;
   double coef_frict = 0.25;
   ContactLaw * contactlaw = new RegularizedCoulombLaw(coef_frict, regularized_time_scale, nb_elements);

/* -------------------------------------------------------------------------- */

   SpectralModel model(nb_elements, nb_time_steps, dom_size, crack_size, nu_mtl, nu_poly, 
		       E_mtl, E_poly, cs_mtl, cs_poly, 
		       tcut_mtl, tcut_poly, overlap, l_index, fracturelaw, contactlaw);
 
   model.initModel();
   model.setLoadingCase(load, load, psi, phi);

   Interfacer interfacer(model);
   interfacer.createCenteredCrack(max_n_str, max_s_str);

   model.updateLoads();
   model.computeInitialVelocities();
   model.computeKernels();

   const std::vector<CrackProfile> & displacements = model.getDisplacements();
   const std::vector<double> & nor_strength = model.getNormalStrength();
   const std::vector<double> & shr_strength = model.getShearStrength();
   const std::vector<Energetics> & E_n = model.getNormalDissipatedEnergy();
   const std::vector<Energetics> & E_s = model.getShearDissipatedEnergy();
   const std::vector<Energetics> & E_fr = model.getFrictionalEnergy();

   int chkpt = 0.75*nb_elements;
   int print = 0.01*nb_time_steps;
   int wdth = 15;
   std::cout  << std::setw(wdth) << "delta_s" 
	      << std::setw(wdth) << "delta_n" 
	      << std::setw(wdth) << "tau_str^s" 
	      << std::setw(wdth) << "tau_str^n" 
	      << std::setw(wdth) << "E_s" 
	      << std::setw(wdth) << "E_n" 
	      << std::setw(wdth) << "E_fr" 
	      << std::endl;

   std::cout.precision(6);

   for (int t = 0; t < nb_time_steps ; ++t) {

     model.updateDisplacements();
     model.updateMaterialProp();
     model.preintegratedKernels();
     model.fftOnDisplacements();
     model.computeTimeConvolution();
     model.computeStresses();
     model.computeVelocities();
     model.computeEnergy();
     model.increaseTimeStep();

     if (print == (int)(0.01*nb_time_steps)) {
    
       std::cout << std::setw(wdth) << displacements[0][chkpt*3]-displacements[1][chkpt*3]  
		 << std::setw(wdth) << displacements[0][chkpt*3+1]-displacements[1][chkpt*3+1] 
		 << std::setw(wdth) << shr_strength[chkpt] 
		 << std::setw(wdth) << nor_strength[chkpt] 
		 << std::setw(wdth) << E_s[0].E 
		 << std::setw(wdth) << E_n[0].E 
		 << std::setw(wdth) << E_fr[0].E 
		 << std::endl;

       print=0;
     }

     ++print;
   }

   delete fracturelaw;
   delete contactlaw;
   
   return 0;
 }

