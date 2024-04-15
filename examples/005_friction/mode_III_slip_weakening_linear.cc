/**
 * @file   mode_III_slip_weakening.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Mar 19 15:47:26 2018
 *
 * @brief  Mode-III crack tip equation of motion
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
#include "cohesive_law_all.hh"
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

int main(int argc, char *argv[])
{

  // Note : Construct the pre-integrated material kernels before running this simulation
  // Use "invert_serial.f" to construct kernel files

  std::cout << "./mode_III_slip_weakening <output_folder_name> <crack_length_ratio> <load_control>" << std::endl;

  std::string output_folder = argv[1];
  // Geometry description

  // UInt nb_time_steps = std::atoi(argv[3]);
  // UInt nex = std::atoi(argv[2]);
  Real cr_ratio = std::atof(argv[2]);
  Real mu = 3e9;
  Real rho = 1200;
  Real nu = 0.33;
  Real E = 2 * mu * (1 + nu);
  Real cs = sqrt(mu / rho);
  Real load_ctrl = std::atof(argv[3]);
  std::cout << "cs = " << cs << std::endl;
  // Cut of the loaded material kernels
  UInt tcut = 100;

  // Loading case
  // Real load = 3e6;
  Real psi = 90.0;
  Real phi = 90.0;

  // Cohesive parameters
  Real crit_n_open = 2e-5;
  Real crit_s_open = 2e-5;
  Real max_n_str = 5e6;
  Real max_s_str = 5e6;
  Real res_n_str = 0.25e6;
  Real res_s_str = 0.25e6;
  Real load = load_ctrl * max_s_str;

  Real G_length = 2 * mu * crit_n_open * (max_n_str - res_n_str) / ((load - res_n_str) * (load - res_n_str) * M_PI);
  // Real G_length = 4*mu*Gc/(M_PI*std::pow(load-res_s_str, 2));

  std::cout << "G_length =" << G_length << std::endl;

  Real dom_sizex = 15 * G_length;

  Real crack_size = cr_ratio * G_length;

  Real lpz = mu * crit_n_open * (max_n_str - res_n_str) / (max_n_str * max_n_str);
  UInt n_ele_ind = std::round(dom_sizex / lpz) * 20;
  UInt nex = 2000;
  UInt nb_time_steps = 5 * nex;
  Real dx = dom_sizex / (Real)(nex);

  std::string sim_name = "Mode-III crack tip equation of motion";

  std::cout << "./mode_III_rate_and_state "
            << "output folder: " << output_folder << "\n"
            << "nb_elements alog x: " << nex << "\n"
            << "nb_time_steps: " << nb_time_steps << "\n"
            << "griffith crack length: " << G_length << "\n"
            << "reference number of elements: " << n_ele_ind << '\n'
            << "crack length ratio: " << cr_ratio << '\n'
            << std::endl;

  /* -------------------------------------------------------------------------- */

  UInt t = 0;
  UInt x_tip = 0;
  UInt x_lap = 0.05 * nex;

  SpectralModel *model;

  model = new SpectralModel(nex, 0, dom_sizex,
                            nu, E, cs, tcut,
                            sim_name, output_folder);

  // Real beta=0.002;
  // SimulationDriver sim_driver(*model, beta=beta);
  SimulationDriver sim_driver(*model);

  Interfacer<_coupled_cohesive> interfacer(*model);

  DataRegister::registerParameter("critical_normal_opening", crit_n_open);
  DataRegister::registerParameter("critical_shear_opening", crit_s_open);
  DataRegister::registerParameter("max_normal_strength", max_n_str);
  DataRegister::registerParameter("max_shear_strength", max_s_str);
  DataRegister::registerParameter("res_shear_strength", res_s_str);
  DataRegister::registerParameter("res_normal_strength", res_n_str);
  interfacer.createUniformInterface();

  interfacer.createThroughCrack((dom_sizex - crack_size) / 2., (dom_sizex + crack_size) / 2.);

  CohesiveLawAll &cohesive_law = dynamic_cast<CohesiveLawAll &>((model->getInterfaceLaw()));

  cohesive_law.preventSurfaceOverlapping(NULL);

  cohesive_law.initRegularFormulation();
  // cohesive_law.initDualFormulation(nor_op_factor, shr_op_factor, nor_str_factor, shr_str_factor);
  // cohesive_law.initExponentialFormulation();
  // cohesive_law.initMultiFormulation(op_list, str_list);

  sim_driver.initConstantLoading(load, psi, phi);

  /* -------------------------------------------------------------------------- */
  // Set-up simulation  outputs

  // DataDumper dumper(*model);
  DataDumper dumper(*model);

/*   dumper.initVectorDumper("ST_Diagram_top_z_velo.cra", _top_velocities, 2, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_z_displ.cra", _top_displacements, 2, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_x_velo.cra", _top_velocities, 1, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_top_x_displ.cra", _top_displacements, 1, 1.0, 1, 0, _text);
  dumper.initVectorDumper("ST_Diagram_shear_stress.cra", _interface_tractions, 2, 1.0, 1, 0, _text);
  dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1.0, 1, 0, _text);
  dumper.initDumper("ST_Diagram_normal_stress.cra", _interface_tractions, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_maximum_shear_strength.cra", _maximum_shear_strength, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_maximum_normal_strength.cra", _maximum_normal_strength, 1.0, 1, 0);
  dumper.initDumper("ST_Diagram_shear_vel.cra", _shear_velocity_jumps, 1.0, 1, 0);
  dumper.initVectorDumper("ST_Diagram_shear_tra.cra", _interface_tractions, 2, 1.0, 1, 0);
  dumper.initVectorDumper("ST_Diagram_normal_tra.cra", _interface_tractions, 1, 1.0, 1, 0); */

  /* -------------------------------------------------------------------------- */

  // sim_driver.launchCrack(dom_sizex/2.,1.75*G_length,0.075,false);

  // DataRegister::restart_dir = "restart_files/";
  // model->restartModel();
  std::ofstream outputFile(output_folder + "ST_cra_tip.cra");
  while ((t < nb_time_steps) && (x_tip < 0.6 * nex))
  {

    // model->pauseModel();
    // break;

    // sim_driver.solveStep();
    model->updateDisplacements();
    model->fftOnDisplacements();
    model->computeStress();
    model->computeInterfaceFields();
    model->increaseTimeStep();

    x_tip = model->getCrackTipPosition(nex / 2, nex);

    if (t % 10 == 0)
    {
      dumper.dumpAll();
      outputFile << x_tip / (Real)(nex) << std::endl;
    }

    if ((x_tip > x_lap) || (t % (UInt)(0.05 * nb_time_steps) == 0))
    {
      std::cout << "Process at " << (Real)t / (Real)nb_time_steps * 100 << "% " << std::endl;
      std::cout << "Crack at " << 100 * x_tip / (Real)(nex) << "% " << std::endl;
      std::cout << std::endl;

      if (x_tip > x_lap)
        x_lap += 0.05 * nex;
    }

    ++t;
  }
  outputFile.close();
  model->pauseModel();

  // delete model;
  return 0;
}
