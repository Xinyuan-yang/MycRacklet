/**
 * @file   cst_speed_fracture.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Nov 11 09:19:07 2015
 *
 * @brief  Dynamic fracture at a constant crack speed in presence of
heterogeneities
 *
 * @section LICENSE
 *
 * cRacklet - A spectral boundary integral method for interface fracture
simulation
 * Copyright (©) 2012 - 2013 Fabian Barras
 *               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 *
 * cRacklet is the result of a collaboration between the Computational Solid
Mechanics
 * Laboratory (LSMS) of Ecole Polytechnique Fédérale de Lausanne (EPFL),
Switzerland
 * and the Department of Aerospace Engineering of the University of Illinois at
 * Urbana-Champaign, United States of America.
 *
 * cRacklet is free software: you can redistribute it and/or modify it under the
terms
 * of the GNU General Public License as published by the Free Software
Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * cRacklet is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
this program.
 * If not, see <http://www.gnu.org/licenses/>.

/* -------------------------------------------------------------------------- */

#include "cohesive_law.hh"
#include "coulomb_law.hh"
#include "data_dumper.hh"
#include "interfacer.hh"
#include "regularized_coulomb_law.hh"
#include "simulation_driver.hh"
#include "spectral_model.hh"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <sys/stat.h>

/* -------------------------------------------------------------------------- */


int main(int argc, char *argv[]) {

  // Note : Construct the pre-integrated material kernels before running this
  // simulation
  // Use "invert_serial.f" to construct kernel files
  // Required command line arguments. [...] delimitates the optional ones

  std::cout << " ./cst_speed_fracture <crack_speed> <loading_file> "
    "<load_writing?> [ <output_folder_name>=./ <nb_t_steps>=10000 "
    "<nb_ele>=8192 <ratio>=0 <loading_angle>=90 ]"
            << std::endl;

  std::string sim_name =
    "Single material homogeneous interface controled rupture speed";

  // Geometry description
  UInt nb_time_steps = 10000;
  UInt nb_elements = 8192;
  Real nu = 0.35;
  Real E = 5.3e9;
  Real cs = 1263;

  // Cut of the loaded material kernels
  UInt tcut = 100;

  // Loading case
  Real load = 4e6;
  Real psi = 90;
  Real phi = 0;

  // Target crack speed v=cr_speed*c_s
  Real cr_speed = std::atof(argv[1]);
  // Name of the loading file
  std::string load_file = argv[2];
  // True=simulation to tailor the loading file
  // False=simulation reading the loading file
  bool write = (bool)(std::atoi(argv[3]));

  // Output folder
  std::string output_folder = "./";
  if (argc > 4)
    output_folder = argv[4];

  if (argc > 5)
    nb_time_steps = std::atoi(argv[5]);

  if (argc > 6)
    nb_elements = std::atoi(argv[6]);

  //ratio asperity strength/surrounding material strength
  Real ratio = 1;

  if (argc > 7)
    ratio = std::atof(argv[7]);

  if (argc > 8)
    psi = std::atof(argv[8]);

  std::cout << "./cst_speed_fracture "
            << "v_imposed:" << cr_speed << " "
            << "load_file:" << load_file << " "
            << "output folder:" << output_folder << " "
            << "nb_time steps:" << nb_time_steps << " "
            << "nb_elements:" << nb_elements << " "
            << "toughness ratio:" << ratio << " "
            << "loading angle:" << psi << " " << std::endl;

  // Cohesive paramters
  Real crit_n_open = 0.01e-3;
  Real crit_s_open = 0.01e-3;
  Real max_s_str = 8e6;
  Real max_n_str = 8e6;
  //Real wk_max_s_str = 4e6;
  //Real wk_max_n_str = 4e6;
  //Real sg_max_n_str = 14e6;
  //Real sg_max_s_str = 14e6;

  Real G_length ;
  if (psi == 45){
    //mixed mode
    G_length =
      sqrt(2) * crit_s_open * max_s_str / (load * load * M_PI) * E / (1 - nu * nu);
  } else {
    //modeII
    G_length =
      crit_s_open * max_s_str / (load * load * M_PI) * E / (1 - nu * nu);
  }


  FractureLaw *fracturelaw;

  std::cout << "Griffith length: " << G_length << std::endl;

  Real dom_sizex = 0.5;
  Real dom_sizez = 0.5 * dom_sizex;
  Real dx = dom_sizex / double(nb_elements);
  UInt  nb_elementsZ = 0.5 * nb_elements ;
  Real dz = dx ;
  Real wall_position = 0.80 * dom_sizex;
  UInt propagation_domain = 0.1 * nb_elements;

  //asperity parameters
  Real asper_radius= 0.02 * dom_sizez ;
  Real asperx = 0.3 * dom_sizex;
  Real asperz = 0.5*dom_sizez;

  // Friction paramters
  bool overlap = 0;
  Real regularized_tcime_scale = 0.1;
  Real coef_frict = 0.25;
  ContactLaw *contactlaw = NULL; //new RegularizedCoulombLaw(
  //coef_frict, regularized_time_scale, nb_elements);

  //phase 0 = 2d initialisation, phase 1 = 3d with asperity
  UInt nb_simulation_phases;
  if(ratio==1)
    nb_simulation_phases=1;
  else
    nb_simulation_phases=2;

  //output_folder of 2d simulation
  std::string outfolder;

  /* --------------------------------------------------------------------------
   */



  for (UInt phase = 0; phase < nb_simulation_phases; ++phase) {

    SpectralModel *model;

    if ((ratio==1.)||(phase==0)){
      outfolder = output_folder + "2d_outputs/";
      mkdir(outfolder.c_str(),0777);

      model = new  SpectralModel({nb_elements, 1}, nb_time_steps, {dom_sizex, 0}, nu, nu,
                                 E, E, cs, cs, tcut, tcut, overlap, fracturelaw,
                                 contactlaw, sim_name, outfolder);
    } else {
      model = new SpectralModel({nb_elements, nb_elementsZ}, nb_time_steps, {dom_sizex, dom_sizez}, nu, nu,
                                E, E, cs, cs, tcut, tcut, overlap, fracturelaw,
                                contactlaw, sim_name, output_folder);
    }







    // SimulationDriver object helping to launch controlled-speed simulation
    SimulationDriver sim_driver(*model, cr_speed, 0.0);

    Interfacer<_linear_coupled_cohesive> interfacer(*model);

    if ((ratio==1.)||(phase==0)){
      interfacer.createUniformInterface(crit_n_open, crit_s_open, max_n_str,
                                        max_s_str);
      interfacer.createThroughCrack(0., 5 * dx);
    } else {
      interfacer.createRightPropagatingCrackRoundAsp (5 * dx, crit_n_open,crit_s_open,max_n_str,
                                                      max_s_str,asper_radius, {asperx,asperz},
                                                      sqrt(ratio),sqrt(ratio));
    }


    interfacer.createThroughWall(wall_position, dom_sizex);
    interfacer.applyInterfaceCreation();

    if (write) {
      // Init algorithm to tailor constant speed loading conditions
      // sim_driver.initConstantSpeed(load, psi, phi, max_s_str);
      sim_driver.initConstantSpeed(load, psi, phi, max_s_str, 0.*dom_sizex, _space_control, 0.9, 0.5*G_length);
    } else {
      // Read an existing file to set loading conditions
      sim_driver.initLoadingFromFile(output_folder + load_file, _space_control,
                                     load, psi, phi);
    }

    DataDumper dumper(*model);

    UInt integ_width = G_length / dx;

    //dumper.initIntegratorsDumper("Energies.cra", {5 * G_length, 1.},
    //                             {10 * G_length, 1.});
    //dumper.initSurfingIntegratorsDumper("surfing_eq.cra", integ_width, 0,
    //                                    nb_elements, {_radiated_energy},
    //                                    {"surfing_Eq"});

    //if (write) {
    dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1, 1, 0, _text);
    // dumper.initDumper("ST_Diagram_shear_velo_jump.cra", _shear_velocity_jumps,
    //                  1, 1, 0, _binary);
    //} else {
    //dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1, 1, 0, _binary);
    //dumper.initDumper("ST_Diagram_shear_velo_jump.cra", _shear_velocity_jumps,
    //                   1, 1, 0, _binary);
    //}

    if(phase==1) {
      DataRegister::restart_dir = "2d_outputs/restart_files/";
      model->restartModel(true);
    }


    if(phase==0) {
      sim_driver.launchCrack(0., G_length, 0.05);
    }

    UInt x_tip = 0;
    //UInt x_tip_prev = x_tip;
    //std::vector<int> iterations(nb_elements);


    UInt t = 0;

    UInt max_t_step = 1.5 * propagation_domain / (model->getBeta() * cr_speed);

    UInt print_info = 0.05 * max_t_step;

    // std::ofstream of;
    // of.open("x_tip.cra");

    while ((x_tip < (propagation_domain)) && (t < max_t_step)) {

      if (phase==0) {
        model->pauseModel();
        std::cout << "End of pseudo 2d" << std::endl;
        break;
      }

      // High level method embedding the resolution of one time step by the
      // SpectralModel
      x_tip = sim_driver.solveStep();

      dumper.dumpAll();

      //if (x_tip == x_tip_prev) {
      //  iterations[x_tip] += 1;
      //} else {
      //  iterations[x_tip] = 1;
      //}

      if (t % ((UInt)print_info) == 0) {
        std::cout << "Crack tip at " << (Real)(x_tip) / (Real)nb_elements * 100
                  << "% "
                  << " or double max time " << (Real)t / (Real)max_t_step * 100
                  << "% " << std::endl;
      }
      //x_tip_prev = x_tip;
      //of << x_tip << std::endl;
      ++t;
    }
    //of.close();

    //
    /*std::ofstream os;
      os.open("iterations.cra");
      for (UInt i = 0; i < iterations.size(); ++i) {
      os << i << " " << iterations[i] << std::endl;
      }
      os.close();*/

    // Do not forget to create the resulting loading file
    // if (write)
    // sim_driver.writeLoading(output_folder + load_file);
    delete model;
  }
  return 0;
}
