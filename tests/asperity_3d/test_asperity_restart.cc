/**
 * @file   test_asperity_restart.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Oct 18 11:11:22 2017
 *
 * @brief  Test asperity after restart
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
/* -------------------------------------------------------------------------- */

 int main(){

   // Note : Construct the pre-integrated material kernels before running this simulation
   // Use "invert_serial.f" to construct kernel files

   // Geometry description
   UInt nb_time_steps =600; 
   UInt nex = 64;
   UInt nez = nex/4;
   Real dom_sizex = 1.0;
   Real dom_sizez = dom_sizex/4;
   Real crack_size = 0.2*dom_sizex;
   Real nu_mtl =  0.33;
   Real nu_poly = 0.35;
   Real E_mtl = 71e9;
   Real E_poly = 53e9;
   Real cs_mtl = 3100;
   Real cs_poly = 1263; 

   // Cut of the loaded material kernels
   UInt tcut_mtl = 100; //mis dans les kernels setup
   UInt tcut_poly = 100;
  
   // Loading case
   Real load = 3e6;
   Real psi = 0;
   Real phi = 0;
   
   // Cohesive parameters
   Real crit_n_open = 0.02e-3;
   Real crit_s_open = 0.02e-3;
   Real max_s_str = 5e6;
   Real max_n_str = 5e6;
   Real ratio = 1.5;

   // Friction parameters
   Real regularized_time_scale = 0.1;
   Real coef_frict = 0.25;

   //output parameters
   UInt direction=0; // if x axis 1 if z axis
   UInt line=0.5*nez;//nex/2-1; //line is a raw or a line and begin with 0
   // Output point position on the line choosen
   std::vector<Real> point_his(8);
   Real dom_size;
   if (direction==0) {dom_size=dom_sizex;}
   else{dom_size=dom_sizez;}
 
   point_his = {dom_size*0.1, dom_size*0.2, dom_size*0.25, dom_size*0.3, dom_size*0.35, dom_sizex*0.4,
		dom_size*0.45, dom_size*0.55, dom_size*0.6, dom_size*0.65, dom_size*0.7, dom_size*0.75, 
		dom_size*0.8, dom_size*0.85, dom_size*0.9};

/* -------------------------------------------------------------------------- */

   UInt print = 0.1*nb_time_steps;
   
   for (UInt step = 0; step < 2; ++step) {
     
     std::shared_ptr<ContactLaw> contactlaw = std::make_shared<RegularizedCoulombLaw>(coef_frict, regularized_time_scale, nex*nez);
     
     SpectralModel model({nex,nez}, nb_time_steps, {dom_sizex,dom_sizez}, 
			 nu_mtl, nu_poly, E_mtl, E_poly, cs_mtl, cs_poly, 
			 tcut_mtl, tcut_poly, "Test asperity");
  
     model.initModel(0.37);
     model.setLoadingCase(load, psi, phi);

     Interfacer<_linear_coupled_cohesive> interfacer(model);
   
     interfacer.createRightPropagatingCrackRoundAsp(crack_size, crit_n_open, crit_s_open, 
						    max_n_str, max_s_str, dom_sizez/8,
						    {0.4*dom_sizex,0.5*dom_sizez}, ratio);
     interfacer.createThroughWall(0.8,1.0);

     std::shared_ptr<CohesiveLaw> cohesive_law = std::dynamic_pointer_cast<CohesiveLaw>(model.getInterfaceLaw());
     cohesive_law->preventSurfaceOverlapping(contactlaw);
     
     model.updateLoads();
     model.initInterfaceFields();
  
     const CrackProfile * t_displacements = model.readData(_top_displacements);
     const CrackProfile * b_displacements = model.readData(_bottom_displacements);
     const std::vector<Real> * nor_strength = model.readData(_normal_strength);
     const std::vector<Real> * shr_strength = model.readData(_shear_strength);

     std::vector<UInt> integ_points_left(nez*nex/2);
  
     UInt local_index = 0;
     for (UInt ix = 0; ix < nex/2; ++ix) {
       for (UInt iz = 0; iz < nez; ++iz) {
	 UInt l_index = ix+iz*nex;
	 integ_points_left[local_index] = l_index;
	 ++local_index;
       }
     }

     std::vector<Real> dx = model.getElementSize();
  
     Integrator E_s(integ_points_left, _shear_fracture_energy, 0., dx[0]*dx[1]);
     model.registerComputer("efrac_shear_left",&E_s);
     Integrator E_n(integ_points_left, _normal_fracture_energy, 0., dx[0]*dx[1]);
     model.registerComputer("efrac_normal_left",&E_n);
     Integrator E_fr(integ_points_left, _frictional_energy, 0., dx[0]*dx[1]);
     model.registerComputer("efric_left",&E_fr);

   
     UInt chkptx = 0.5*nex;
     UInt chkptz = 0.5*nez;
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

     UInt t_start, t_end;
          
     if(step==1) {
       model.restartModel();
       t_start = 301;
       t_end = nb_time_steps;
     }
     else {
       t_start = 0;
       t_end = 301;
     }
       	 
     for (UInt t = t_start; t < t_end; ++t) {
     
       model.updateDisplacements();
       model.fftOnDisplacements();
       model.computeStress();
       model.computeInterfaceFields();
       model.increaseTimeStep();

       //if(t==301)
       //
     
       if (print == (UInt)(0.1*nb_time_steps)) {
	 std::cout << std::setw(wdth) << (*t_displacements)[(chkptx-1+nex*(chkptz-1))*3+2]-(*b_displacements)[(chkptx-1+nex*(chkptz-1)
													       )*3+2]  
		   << std::setw(wdth) << (*t_displacements)[(chkptx-1+nex*(chkptz-1))*3+1]-(*b_displacements)[(chkptx-1+nex*(chkptz-1
															     ))*3+1] 
		   << std::setw(wdth) << (*shr_strength)[(chkptx-1+nex*(chkptz-1))] 
		   << std::setw(wdth) << (*nor_strength)[(chkptx-1+nex*(chkptz-1))] 
		   << std::setw(wdth) << E_s.getIntegration()  
		   << std::setw(wdth) << E_n.getIntegration() 
		   << std::setw(wdth) << E_fr.getIntegration() 
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
