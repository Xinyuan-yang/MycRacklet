/**
 * @file   alu_homa3d.cc
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Nov 16 10:03:54 2015
 *
 * @brief  3d interface in the presence of one asperity (See http://lsms.epfl.ch/page-74218-en.html)
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
/* -------------------------------------------------------------------------- */

 int main(int argc, char *argv[]){

   // Note : Construct the pre-integrated material kernels before running this simulation
   // Use "invert_serial.f" to construct kernel files

  std::cout << "./single_asp_3d <output_folder_name> <asperity_toughness ratio> <nb_ele_x> <nb_time_steps> [ <loading_angle=0.0> ]" << std::endl;

   std::string output_folder=argv[1];
   
   // Geometry description
   UInt nb_time_steps = std::atoi(argv[4]); 
   UInt nex = std::atoi(argv[3]); 
   UInt nez = nex/4;
   Real nu =  0.35;
   Real E = 5.3e9;
   Real cs = 1263;

   // Cut of the loaded material kernels
   UInt tcut = 100; 
  
   // Loading case
   Real load = 1e6;
   Real psi =0.0;
   if(argc > 5)
     psi = std::atof(argv[5]);
   Real phi = 0.0;

   // Cohesive parameters
   Real crit_n_open = 0.02e-3;//0.08e-3;
   Real crit_s_open = 0.02e-3;//0.08e-3;
   Real max_n_str = 5e6;//1.25e6;
   Real max_s_str = 5e6;//1.25e6;
   Real G_length = crit_n_open*max_n_str/(load*load*M_PI)*E/(1-nu*nu);

   Real dom_sizex = 24*G_length;
   Real dom_sizez = 3*G_length;
   Real dx = dom_sizex/(Real)(nex);
   Real dz = dom_sizez/(Real)(nez);

   Real crack_size = 5*dx;//0.5*G_length;
   
   //ratio asperity strength/surrounding material strength
   Real ratio = 5;
   ratio = std::atof(argv[2]);
   Real asper_radius=0.5*G_length;
   Real asperx = 2*G_length;
   Real asperz = 0.5*dom_sizez;
   FractureLaw * fracturelaw;

   // Friction parameters
   bool overlap = 0;
   Real regularized_time_scale = 0.1;
   Real coef_frict = 0.25;
   ContactLaw * contactlaw = new RegularizedCoulombLaw(coef_frict, regularized_time_scale, nex*nez);

   std::string sim_name = "Mode-I fracture along a 3d interface with a single asperity";

   std::cout << "./single_asp_3d " 
	     << "output folder: " << output_folder << " " 
	     << "toughness ratio: " << ratio << " " 
	     << "nb_elements alog x: " << nex << " "
	     << "nb_time_steps: " << nb_time_steps << " "
	     << "griffith crack length: " << G_length << " " 
	     << std::endl;

   
/* -------------------------------------------------------------------------- */
   SpectralModel * model;

   if (ratio==1.)
     model = new SpectralModel({nex,1}, nb_time_steps, {dom_sizex,0.},
			       nu, nu, E, E, cs, cs, tcut, tcut, overlap,
			       fracturelaw, contactlaw, sim_name, output_folder);
   else
     model = new SpectralModel({nex,nez}, nb_time_steps, {dom_sizex,dom_sizez},
			       nu, nu, E, E, cs, cs, tcut, tcut, overlap,
			       fracturelaw, contactlaw, sim_name, output_folder);

   SimulationDriver sim_driver(*model);
   
   Interfacer<_linear_coupled_cohesive> interfacer(*model);   

   if (ratio==1.) {
     interfacer.createUniformInterface(crit_n_open, crit_s_open, 
				       max_n_str, max_s_str);
     interfacer.createThroughCrack(0.,crack_size);
   } else
     interfacer.createRightPropagatingCrackRoundAsp(crack_size, crit_n_open,crit_s_open,max_n_str,
						    max_s_str,asper_radius, {asperx,asperz},
						    sqrt(ratio),sqrt(ratio));

   interfacer.createThroughWall(0.95*dom_sizex,dom_sizex);
   interfacer.applyInterfaceCreation();
   
   sim_driver.initConstantLoading(load, psi, phi);

   DataDumper dumper(*model);

   /* -------------------------------------------------------------------------- */

   std::vector<Real> start_corner = {0.4,0.15};
   std::vector<Real> end_corner = {1.0,0.35};
   
   //dumper.initIntegratorsDumper("Energy_local.cra",start_corner,end_corner);
   
   /* -------------------------------------------------------------------------- */

   Real x_start = 0.4;
   Real x_end = 1.0;
   UInt step = 5;
   std::vector<Real> z_coord;
   
   if (ratio==1.)
     z_coord = {0.};
   else
     z_coord = {0.25, 0.3, 0.35, 0.4};
   
   std::vector<UInt> obs_points;

   for (UInt z = 0; z < z_coord.size(); ++z) {
     UInt iz = (z_coord[z]/dz);
     UInt ix = (x_start/dx);
     while (dx*(Real)ix<x_end){
       obs_points.push_back(iz*nex+ix);
       ix += step;
     }
   }
   std::vector<DataFields> fields =
     {_interface_tractions,_normal_strength,_normal_displacement_jumps,_normal_velocity_jumps};

   std::cout << obs_points.size() << " observation points used to monitor fields evolution." << std::endl;   
   //dumper.initPointsDumper("mid-points_fields.cra", fields, obs_points, _binary);

   /* -------------------------------------------------------------------------- */

   //dumper.initIntegratorsDumper("Energy_global.cra");
   dumper.initDumper("ST_Diagram_nor_velo_jump.cra", _normal_velocity_jumps, 0.25, 4, 0, _binary);
   dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1.0, 1, 0, _binary);
   dumper.initDumper("ST_Diagram_nor_strength.cra", _normal_strength, 0.25, 4, 0, _binary);
  
   sim_driver.launchCrack(0.,2*G_length,0.05);
   
   UInt x_tip=0;
   UInt x_lap = 0.05*nex;
   for (UInt t = 0; t < nb_time_steps ; ++t) {

     sim_driver.solveStep();
     x_tip = model->getCrackTipPosition(0.,nex);
     
     if (t%100==0)
       dumper.dumpAll();
     //       dumper.dump("mid-points_fields.cra");
       //dumper.dump("Energy_local.cra");
          
     if ((x_tip>x_lap)||(t%(UInt)(0.05*nb_time_steps)==0)) {
       std::cout << "Process at " << (Real)t/(Real)nb_time_steps*100 << "% " << std::endl;
       std::cout << "Crack at " << 100*x_tip/(Real)(nex) << "% " << std::endl;
       if (x_tip>x_lap)
	 x_lap += 0.05*nex;
     }
   }

   delete contactlaw;
   
   return 0;
 }

