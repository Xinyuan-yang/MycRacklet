/**
 * @file   cst_speed_fracture_2d_3d.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 * @author Fatima Fekak <fatima-ezzahra.fekak@epfl.ch>
 * @date   Wed Nov 11 09:19:07 2015
 *
 * @brief  Dynamic fracture at a constant crack speed in presence of heterogeneities
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
#include "data_dumper.hh"
#include "interfacer.hh"
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

  std::cout << " ./cst_speed_fracture_2d_3d <bool=is it a load calibrating simulation?> "
    "<input_file_name> <output_folder_name>"
    " <bool; set 1 to run a prediction on memory requirement (default=0)>" << std::endl;

  std::string sim_name =
    "Crack front waves along a constant speed rutpure front";

  // True= 2D calibration simulation to tailor the loading file
  // False= 3D simulation reading the loading file
  bool write = (bool)(std::atoi(argv[1]));
   
  std::string input_file = argv[2];
  DataRegister::readInputFile(input_file);

  // Output folder
  std::string output_folder = argv[3];

  bool memory_test = false;
  
  if (argc > 4)
    memory_test = argv[4];
  
  // Geometry description
  UInt nb_elements = DataRegister::getParameter<UInt>("nb_elements");
  Real nu = DataRegister::getParameter<Real>("nu");
  Real E = DataRegister::getParameter<Real>("E");
  Real cs = DataRegister::getParameter<Real>("c_s");

  // Loading case
  Real load = DataRegister::getParameter<Real>("far_field_load");
  Real psi = DataRegister::getParameter<Real>("load_angle_psi");
  Real phi = DataRegister::getParameter<Real>("load_angle_phi");

  // Target crack speed v=cr_speed*c_s
  Real cr_speed = DataRegister::getParameter<Real>("constant_crack_speed");

  //ratio asperity strength/surrounding material strength
  Real ratio = DataRegister::getParameter<Real>("asperity_toughness_ratio");
  
  // Cohesive paramters
  Real crit_n_open  = DataRegister::getParameter<Real>("critical_normal_opening");
  Real crit_s_open = DataRegister::getParameter<Real>("critical_shear_opening");
  Real max_s_str = DataRegister::getParameter<Real>("max_shear_strength");
  Real max_n_str = DataRegister::getParameter<Real>("max_normal_strength");
    
  // Associated griffith length
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
  
  Real dom_sizex = DataRegister::getParameter<Real>("dom_size_x");
  Real dom_sizez = DataRegister::getParameter<Real>("dom_size_z");
  Real dx = dom_sizex / double(nb_elements);
  UInt  nb_elementsZ = nb_elements * dom_sizez/dom_sizex;
  Real dz = dx ;
  UInt propagation_domain = DataRegister::getParameter<Real>("max_crack_size") *nb_elements;
  
  //asperity parameters
  Real asper_radius= DataRegister::getParameter<Real>("asperity_radius") * dom_sizez ;
  Real asperx = DataRegister::getParameter<Real>("asperity_pos_x") * dom_sizex;
  Real asperz = DataRegister::getParameter<Real>("asperity_pos_z") *dom_sizez;

  // Cut off time of the elastodynamic kernels
  UInt tcut = 100;
  
  std::cout << "Brief summary: "
	    << "Calibration simulation ? " << write << std::endl
            << "v_imposed:" << cr_speed << std::endl
            << "output folder:" << output_folder << std::endl
            << "nb_elements:" << nb_elements << std::endl
            << "loading psi angle:" << psi << std::endl
	    << "Domain size: {" << dom_sizex << "," << dom_sizez << "}" << std::endl
	    << "Griffith length: " << G_length << std::endl
	    << "Delta G_c / G_c: " << ratio << std::endl
	    << "Asperity radius: " << asper_radius << std::endl
	    << "Asperity position: {" <<asperx << "," << asperz << "}" << std::endl;
  
  //phase 0 = 2d initialisation, phase 1 = 3d with asperity
  UInt nb_simulation_phases;
  if(write)
    nb_simulation_phases=1;
  else
    nb_simulation_phases=2;

  //output_folder of 2d simulation
  std::string outfolder;
  // Name of the loading file
  std::string load_file = "loading_file.cra";

/* -------------------------------------------------------------------------- */

  for (UInt phase = 0; phase < nb_simulation_phases; ++phase) {

    SpectralModel *model;
    
    // When the number of time steps is not known a priori, setting it to zero will initialize the
    // largest model possible (i.e. allocate enough memory to store displacements history up to the cut-off time)
    UInt nb_time_steps = 0;

    if ((write)||(phase==0)){
      // Create 2d model and the output folder where the loading file will be written
      outfolder = output_folder + "2d_outputs/";
      mkdir(outfolder.c_str(),0777);

      model = new  SpectralModel(nb_elements, nb_time_steps, dom_sizex, nu, 
				 E, cs, tcut, sim_name, outfolder);
    }
    else {
      // Create 3d model
      model = new SpectralModel({nb_elements, nb_elementsZ}, nb_time_steps, {dom_sizex, dom_sizez}, nu,
                                E, cs, tcut, sim_name, output_folder);
    }

    // Blank model initialiasation to predict the required memory size
    if(memory_test) {
      UInt nb_it=std::ceil(1/(cr_speed*CS_DT_OVER_DX));
      Real new_beta = 1.0/(nb_it*cr_speed);    
      model->initModel(new_beta,true);
      continue;
    }
    
    // SimulationDriver helps with running controlled-speed simulation
    SimulationDriver sim_driver(*model, cr_speed, 0.0);

    // Interfacer helps with generating the interface properties
    Interfacer<_linear_coupled_cohesive> interfacer(*model);

    if ((write)||(phase==0)){
      // In the 2d case, uniform interface properties are defined.
      // The values of the uniform properties are directly searched in the input file
      interfacer.createUniformInterface();
      interfacer.createThroughCrack(0., 5 * dx);
    } else {
      // In case of a 3d simulation, a circular asperity stands in the middle of a homogeneous interface
      interfacer.createRightPropagatingCrackRoundAsp (5 * dx, crit_n_open,crit_s_open,max_n_str,
                                                      max_s_str,asper_radius, {asperx,asperz},
                                                      sqrt(ratio),sqrt(ratio));
    }
    // A tougher wall is added at the end of the domain to ensure crack propagation in a single direction
    Real wall_position = 0.80 * dom_sizex;
    interfacer.createThroughWall(wall_position, dom_sizex);

    if (write) {
      // Init algorithm to tailor constant speed loading conditions
      // See simulation_driver.hh for more information
      sim_driver.initConstantSpeed(load, psi, phi, max_s_str, 0.*dom_sizex, _space_control, 0.9, 0.5*G_length);
    } else {
      // Read an existing file to set loading conditions
      sim_driver.initLoadingFromFile(output_folder + load_file, _space_control, load, psi, phi);
    }

    // DataDumper helps with generating output of the simulation    
    DataDumper dumper(*model);

    // Initialize an output of crack index (i.e. 0=intact, 1=cohesive_zone, 2=broken)
    dumper.initDumper("ST_Diagram_id.cra", _id_crack);
    // Check DataFields in data_register.hh to output other kinds of simulation data
    // Check DataDumper in data_dumper.hh to generate different output formats
    // dumper.initDumper("ST_Diagram_shear_velo_jump.cra", _shear_velocity_jumps);

    if(phase==0) {
      // This function slowly grows a seed crack from a given position up to a given crack size
      // This artifical growth occurs at a given rate (e.g. below 0.025*c_s)
      sim_driver.launchCrack(0., G_length, 0.025);
    }

    if(phase==1) {
      // Restart the 3d model from the state obtained after the artifical growth of a seed crack computed in 2d
      DataRegister::restart_dir = "2d_outputs/restart_files/";
      model->restartModel(true);
    }

/* -------------------------------------------------------------------------- */
    UInt x_tip = 0;
    UInt t = 0;
    UInt max_t_step = 1.5 * propagation_domain / (model->getBeta() * cr_speed);
    UInt print_info = 0.05 * max_t_step;
    
    while ((x_tip < (propagation_domain)) && (t < max_t_step)) {

      if ((phase==0)&&(!write)) {
	  // Stop the initialization phase after the artifical growth of a 2d seed crack and the simulation state
	  model->pauseModel();
	  std::cout << "End of pseudo 2d" << std::endl;
	  break;
      }

      // Method solving one time step by SpectralModel and returning the current position of the crack tip
      x_tip = sim_driver.solveStep();

      // The current state of simulation fields is printed for every initialized output files
      dumper.dumpAll();

      if (t % ((UInt)print_info) == 0) {
        std::cout << "Crack tip at " << (Real)(x_tip) / (Real)nb_elements * 100
                  << "% "
                  << " or double max time " << (Real)t / (Real)max_t_step * 100
                  << "% " << std::endl;
      }
      ++t;
    }
    
    // The loading file is written at the end of the 2d calibration of a constant speed fracture
    if (write)
      sim_driver.writeLoading(output_folder + load_file);
    
    delete model;
  }
  
  return 0;
}
