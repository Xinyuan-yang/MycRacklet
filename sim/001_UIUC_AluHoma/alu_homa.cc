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
/* -------------------------------------------------------------------------- */

 int main(){

   // Note : Construct the pre-integrated material kernels before running this simulation
   // Use "invert_serial.f" to construct kernel files

   // Geometry description
   int nb_time_steps = 7000; 
   int nb_elements = 4096;
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

   // Output point position
   std::vector<double> point_his(8);
   point_his = {dom_size*0.1, dom_size*0.2, dom_size*0.25, dom_size*0.3, dom_size*0.35, dom_size*0.4,
		dom_size*0.45, dom_size*0.55, dom_size*0.6, dom_size*0.65, dom_size*0.7, dom_size*0.75, 
		dom_size*0.8, dom_size*0.85, dom_size*0.9};

/* -------------------------------------------------------------------------- */

   SpectralModel model(nb_elements, nb_time_steps, dom_size, crack_size, nu_mtl, 
		       nu_poly, E_mtl, E_poly, cs_mtl, cs_poly, 
		       tcut_mtl, tcut_poly, overlap, l_index, fracturelaw, contactlaw);
  
   model.initModel();
   model.setLoadingCase(load, load, psi, phi);

   Interfacer interfacer(model);
   interfacer.createCenteredCrack(max_n_str, max_s_str);

   model.updateLoads();
   model.computeInitialVelocities();
   model.computeKernels();

   extern std::string lm_release_info;

   DataDumper dumper(model, "Mixed-mode debonding at Aluminium-Homalite interface", lm_release_info);
   dumper.initEnergetics("Energy.dat");
   dumper.initSpaceTimeDiagram("ST_Diagram_id.dat", _cracking_index);
   dumper.initSpaceTimeDiagram("ST_Diagram_normal.dat", _normal_traction);
   dumper.initSpaceTimeDiagram("ST_Diagram_shear.dat", _shear_traction);
   dumper.initSnapshot("top_displ_snapshot.dat", _top_displacement);
   dumper.initSnapshot("bot_displ_snapshot.dat", _bot_displacement);
   dumper.initSnapshot("tract_snapshot.dat", _tractions);
   dumper.initPointsHistory("Points_history.dat", point_his);

   int print = 0.1*nb_time_steps;

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

     dumper.printEnergetics(t);
     dumper.printSpaceTimeDiagram(t);
     dumper.printPointsHistory(t);

     if (print == (int)(0.1*nb_time_steps)) {
       std::cout << "Process at " << (double)t/(double)nb_time_steps*100 << "% " << std::endl;
       print=0;

       dumper.printSnapshot(t);
     }

     ++print;
   }


   delete fracturelaw;
   delete contactlaw;
   
   return 0;
 }

