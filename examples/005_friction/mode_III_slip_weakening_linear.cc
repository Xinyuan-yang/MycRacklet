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

  std::cout << "./mode_III_slip_weakening <output_folder_name> <nb_ele_x> <nb_time_steps> <crack_length_ratio>" << std::endl;

  std::string output_folder = argv[1];
  // Geometry description

  UInt nb_time_steps = std::atoi(argv[3]);
  UInt nex = std::atoi(argv[2]);
  Real cr_ratio = std::atof(argv[4]);
  Real mu = 3e9;
  Real rho = 1200;
  Real nu = 0.33;
  Real E = 2 * mu * (1 + nu);
  Real cs = sqrt(mu / rho);
  std::cout << "cs = " << cs << std::endl;
  bool incre = true;
  // Cut of the loaded material kernels
  UInt tcut = 100;

  // Loading case
  Real load = 0.8e6;
  Real psi = 90.0;
  Real phi = 90.0;

  // Cohesive parameters
  Real crit_n_open = 50.0e-5;
  Real crit_s_open = 50.0e-5;
  Real max_n_str = 5e6;
  Real max_s_str = 5e6;
  Real res_n_str = 0.00e6;
  Real res_s_str = 0.00e6;

  Real incr_x;
  Real incr_y;
  Real incr_z;
  Real load_actu = 0;

  Real Gc = 0.5*crit_n_open*(max_n_str-res_n_str);
  Real G_length = 4.0 * mu * Gc / (M_PI*std::pow(load-res_s_str, 2));
  //Real G_length = 2*E/(1-nu*nu)*Gc/(M_PI*std::pow(load-res_s_str, 2));

  std::cout << "G_length =" << G_length << std::endl;

  Real dom_sizex = 15 * G_length;
  //dom_sizex = 8;

  Real crack_size = cr_ratio * G_length;

  Real lpz = mu * crit_n_open * (max_n_str - res_n_str) / (max_n_str * max_n_str);
  UInt n_ele_ind = std::round(dom_sizex / lpz) * 20;
  Real dx = dom_sizex / (Real)(nex);

  if(lpz >= G_length/20) std::cerr << "\033[33mWarning : process zone size not small enough compared to the Griffith crack length, Griffith's theory may not be valid\033[0m"<<std::endl;

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

  model = new SpectralModel(nex, nb_time_steps, dom_sizex,
                            nu, E, cs, tcut,
                            sim_name, output_folder);

  // Real beta=0.002;
  // SimulationDriver sim_driver(*model, beta=beta);
  model->initModel();

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

  model->setLoadingCase(load, psi, phi);
  if (incre)
  {
    psi *= M_PI/180;
    phi *= M_PI/180;
    model->setLoadingCase(0, psi, phi);
    incr_x = load * sin(psi) * cos(phi) / nb_time_steps * 2.0;
    incr_y = load * cos(psi) / nb_time_steps * 2.0;
    incr_z = load * sin(psi) * sin(phi) / nb_time_steps * 2.0;
  }
  std::cout<< sin(psi) << std::endl;
  model->updateLoads();
  model->initInterfaceFields();
  x_tip = model->getCrackTipPosition(nex / 2, nex);
  x_lap = x_tip + 0.01 * nex;

  /* -------------------------------------------------------------------------- */
  // Set-up simulation  outputs

  // DataDumper dumper(*model);
  DataDumper dumper(*model);

    dumper.initVectorDumper("ST_Diagram_top_z_velo.cra", _top_velocities, 2, 1.0, 1, 0, _text);
    dumper.initVectorDumper("ST_Diagram_top_z_displ.cra", _top_displacements, 2, 1.0, 1, 0, _text);
    dumper.initVectorDumper("ST_Diagram_top_x_velo.cra", _top_velocities, 1, 1.0, 1, 0, _text);
    dumper.initVectorDumper("ST_Diagram_top_x_displ.cra", _top_displacements, 1, 1.0, 1, 0, _text);
    dumper.initVectorDumper("ST_Diagram_shear_stress.cra", _interface_tractions, 2, 1.0, 1, 0, _text);
    dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1.0, 1, 0, _text);
  // sim_driver.launchCrack(dom_sizex/2.,1.75*G_length,0.075,false);

  // DataRegister::restart_dir = "restart_files/";
  // model->restartModel();
  std::ofstream outputtip(output_folder + "ST_cra_tip.cra");
  std::ofstream outputload(output_folder + "ST_load.cra");
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
    if (incre && (t < nb_time_steps / 2.))
    {
      model->incrementLoad(incr_x, 0);
      model->incrementLoad(incr_y, 1);
      model->incrementLoad(incr_z, 2);
      load_actu += incr_z;
    }

    if (t % 10 == 0)
    {
      dumper.dumpAll();
      outputload <<  load_actu << std::endl;
      outputtip << x_tip / (Real)(nex) << std::endl;
    }

    if ((x_tip > x_lap) || (t % (UInt)(0.05 * nb_time_steps) == 0))
    {
      std::cout << "Process at " << (Real)t / (Real)nb_time_steps * 100 << "% " << std::endl;
      std::cout << "Crack at " << 100 * x_tip / (Real)(nex) << "% " << std::endl;
      std::cout << std::endl;

      if (x_tip > x_lap)
        {x_lap += 0.05 * nex;
      }
    }

    ++t;
  }
  outputtip.close();
  outputload.close();
  model->pauseModel();

  // delete model;
  return 0;
}
