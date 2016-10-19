/**
 * @file   alu_homa3d.cc
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
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

  std::cout << "./single_asp_3d <output_folder_name> <asperity_toughness ratio> <nb_ele_x> <nb_time_steps>" << std::endl;


   std::string output_folder=argv[1];
   
   // Geometry description
   UInt nb_time_steps = std::atoi(argv[4]); 
   UInt nex = std::atoi(argv[3]); 
   UInt nez = nex/2;
   Real dom_sizex = 4.0;
   Real dom_sizez = 2.0;
   Real crack_size = dom_sizex/4;
   Real nu =  0.35;
   Real E = 81e9;
   Real cs = 3000;

   // Cut of the loaded material kernels
   UInt tcut = 100; 
  
   // Loading case
   Real load = 3e6;
   Real psi =89.9;
   Real phi = 0.0;
   UInt l_index = 1;

   // Cohesive parameters
   Real crit_n_open = 0.02e-3;
   Real crit_s_open = 0.02e-3;
   Real max_n_str = 5e6;
   Real max_s_str = 5e6;

   //ratio asperity strength/surrounding material strength
   Real ratio = 5;
   ratio = std::atoi(argv[2]);
   Real asper_radius=0.2;
   Real asperx = 0.4*dom_sizex;
   Real asperz = 0.5*dom_sizez;
   FractureLaw * fracturelaw;

   // Friction parameters
   bool overlap = 0;
   Real regularized_time_scale = 0.1;
   Real coef_frict = 0.25;
   ContactLaw * contactlaw = new RegularizedCoulombLaw(coef_frict, regularized_time_scale, nex*nez);

   std::string sim_name = "Slip along a 3d interface with a single asperity";
   
/* -------------------------------------------------------------------------- */

   SpectralModel model({nex,nez}, nb_time_steps, {dom_sizex,dom_sizez},
		       nu, nu, E, E, cs, cs, tcut, tcut, overlap, l_index,
		       fracturelaw, contactlaw, sim_name, output_folder);

   SimulationDriver sim_driver(model);
   
   Interfacer<_linear_coupled_cohesive> interfacer(model);   
   interfacer.createRightPropagatingCrackRoundAsp(crack_size, crit_n_open,crit_s_open,max_n_str,
						  max_s_str,asper_radius, {asperx,asperz}, ratio);
   interfacer.createThroughWall(0.95*dom_sizex,dom_sizex);
   interfacer.applyInterfaceCreation();
   
   sim_driver.initConstantLoading(load, psi, phi);

   UInt chkpt = 0.25*nex*nex+0.5*nez;
   UInt wdth = 15;

   const CrackProfile * t_displacements = model.readData(_top_displacements);
   const CrackProfile * b_displacements = model.readData(_bottom_displacements);
   DataDumper dumper(model);

   dumper.initEnergetics("Energy3d.cra");
   dumper.initDumper("ST_Diagram_shr_velo_jump.cra", _shear_velocity_jumps);
   dumper.initDumper("ST_Diagram_id.cra", _id_crack);
   dumper.initDumper("ST_Diagram_shr_strength.cra", _shear_strength);
  
   for (UInt t = 0; t < nb_time_steps ; ++t) {

     sim_driver.solveStep();
          
     dumper.printEnergetics();

     if (t%25==0){
     dumper.dumpAll();
     }

     if (t%(UInt)(0.05*nb_time_steps)==0) {
       
       std::cout << "Process at " << (Real)t/(Real)nb_time_steps*100 << "% " << std::endl;
       std::cout << std::setw(wdth) <<  (*t_displacements)[chkpt*3]-(*b_displacements)[chkpt*3] << std::endl; 
     }
   }

   delete contactlaw;
   
   return 0;
 }

