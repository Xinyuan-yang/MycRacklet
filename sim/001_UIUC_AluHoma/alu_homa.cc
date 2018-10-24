/**
 * @file   alu_homa.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Mar 11 11:19:07 2013
 *
 * @brief  Study dynamic debonding at a planar interface btw Aluminium (mtl) and Homalite (poly) 
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
 */
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
   UInt nb_time_steps = 7000; 
   UInt nb_elements = 4096;
   Real dom_size = 1.0;
   Real crack_size = 0.05;
   Real nu_mtl = 0.33;
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
   Real regularized_time_scale = 0.1;
   Real coef_frict = 0.25;
   ContactLaw * contactlaw = new RegularizedCoulombLaw(coef_frict, regularized_time_scale, nb_elements);

   // Output point position
   std::vector<UInt> points_int;
   points_int = {(UInt)(nb_elements*0.1), (UInt)(nb_elements*0.2), (UInt)(nb_elements*0.25), (UInt)(nb_elements*0.3),
		 (UInt)(nb_elements*0.35), (UInt)(nb_elements*0.4), (UInt)(nb_elements*0.45), (UInt)(nb_elements*0.55),
		 (UInt)(nb_elements*0.6), (UInt)(nb_elements*0.65), (UInt)(nb_elements*0.7), (UInt)(nb_elements*0.75),
		 (UInt)(nb_elements*0.8), (UInt)(nb_elements*0.85), (UInt)(nb_elements*0.9) };

   /* -------------------------------------------------------------------------- */
   SpectralModel model({nb_elements,1}, nb_time_steps, {dom_size,0.}, nu_mtl, 
		       nu_poly, E_mtl, E_poly, cs_mtl, cs_poly, 
		       tcut_mtl, tcut_poly, 
		       "Mixed-mode debonding at Aluminium-Homalite interface");

   Real beta=0.4; //Stable time step coefficient
   
   model.initModel(beta);
   model.setLoadingCase(load, psi, phi);

   Interfacer<_linear_coupled_cohesive> interfacer(model);
   interfacer.createThroughCenteredCrack(crack_size, crit_n_open, crit_s_open, max_n_str, max_s_str);

   CohesiveLaw * cohesive_law = dynamic_cast<CohesiveLaw*>(*(model.getInterfaceLaw()));
   cohesive_law->preventSurfaceOverlapping(contactlaw);

   model.updateLoads();
   model.initInterfaceFields();

   DataDumper dumper(model);
   std::string st_diag_id = "ST_Diagram_id.cra";
   std::string st_diag_nor_trac = "ST_Diagram_normal_tractions.cra";
   std::string st_diag_shear_velo = "ST_Diagram_shear_velocity_jumps.cra";
   std::string top_u = "top_displ_snapshot.cra";
   std::string bot_u = "bot_displ_snapshot.cra";
   std::string tractions = "trac_snapshots.cra"; 
   std::string point_his = "Points_history.cra";
   std::string energy = "Energy.cra";
   dumper.initDumper(st_diag_shear_velo, _shear_velocity_jumps);
   dumper.initVectorDumper(st_diag_nor_trac, _interface_tractions,1);
   dumper.initDumper(st_diag_id, _id_crack);   
   dumper.initDumper(top_u, _top_displacements);
   dumper.initDumper(bot_u, _bottom_displacements);
   dumper.initDumper(tractions, _interface_tractions);
   dumper.initPointsDumper(point_his, points_int); 
   dumper.initIntegratorsDumper(energy);
   UInt print = 0.1*nb_time_steps;

   for (UInt t = 0; t < nb_time_steps ; ++t) {

     model.updateDisplacements();
     model.fftOnDisplacements();
     model.computeStress();
     model.computeInterfaceFields();
     model.increaseTimeStep();

     dumper.dump(st_diag_id);
     dumper.dump(st_diag_nor_trac);
     dumper.dump(st_diag_shear_velo);
     dumper.dump(point_his);
     dumper.dump(energy);
     
     if (print == (UInt)(0.1*nb_time_steps)) {
       std::cout << "Process at " << (Real)t/(Real)nb_time_steps*100 << "% " << std::endl;
       print=0;

       dumper.dump(top_u);
       dumper.dump(bot_u);
       dumper.dump(tractions);
     }

     ++print;
   }
   
   return 0;
 }

