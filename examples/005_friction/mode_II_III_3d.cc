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
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>
/* -------------------------------------------------------------------------- */
// fluid overpressure
inline std::vector<Real> fluidoverpressure(std::vector<Real> load, UInt nex, UInt nez, Real delta_p_star, Real alpha, Real dx, Real t)
{
    // To avoid underflow report errors
    Real x_mid = nex / 2.0;
    Real z_mid = nez / 2.0;
    Real r;
    gsl_set_error_handler_off();
    for (UInt x = 0; x < nex; x++)
    {
        for (UInt z = 0; z < nez; z++)
        {
            r = sqrt((x - x_mid) * (x - x_mid) + (z-z_mid)*(z-z_mid) ) * dx;
            if (r == 0)
            {
                r = dx;
            }
            load[x * 3 + z * nex * 3 + 1] += 1 * delta_p_star * gsl_sf_expint_E1(r * r / 4 / alpha / t);
        }
    }

    return load;
}

int main(int argc, char *argv[])
{

    // Note : Construct the pre-integrated material kernels before running this simulation
    // Use "invert_serial.f" to construct kernel files

    std::cout << "./mode_III_slip_weakening <output_folder_name> <nb_ele_x> <nb_time_steps> <bool_exp>" << std::endl;

    std::string output_folder = argv[1];
    // Geometry description

    UInt nb_time_steps = std::atoi(argv[3]);
    UInt nex = std::atoi(argv[2]);
    UInt nez = nex;
    Real mu = 11.84e9;
    Real rho = 1200;
    Real nu = 0.01;
    Real E = 2 * mu * (1 + nu);
    Real cs = sqrt(mu / rho);
    std::cout << "cs = " << cs << std::endl;
    // Cut of the loaded material kernels
    UInt tcut = 100;
    bool exp = std::atoi(argv[4]);

    // Loading case
    Real load_nor = 5.08e6;
    Real load_shr = 1.83e6;
    // Cohesive parameters
    Real crit_n_open = 0.37e-3;
    Real crit_s_open = 0.37e-3;
    Real max_n_str = 5e6;
    Real max_s_str = 5e6;
    Real res_n_str = 0.25e6;
    Real res_s_str = 0.25e6;
    Real delta_p_star = 0.035 * load_nor;
    Real alpha = 0.88e6;
    Real mus = 0.6;
    Real mud = 0.42;

    Real G_length = 2 * mu * crit_n_open * load_nor * (mus - mud) / (load_shr * load_shr * (1 - mud) * (1 - mud) * M_PI);
    Real R_w = mu * crit_n_open / (mus - mud) / load_nor;
    // Real G_length = 4*mu*Gc/(M_PI*std::pow(load-res_s_str, 2));

    Real dom_sizex = 20 * R_w;
    Real dom_sizez = dom_sizex;
    Real dx = dom_sizex / (Real)(nex);
    Real crack_size = 0 * G_length;
    // Real crack_size = 0.1;

    Real lpz = mu * crit_n_open * (max_s_str - res_s_str) / (max_s_str * max_s_str);
    UInt n_ele_ind = std::round(dom_sizex / lpz) * 20;

    std::string sim_name = "Mode-III crack tip equation of motion";

    std::cout << "./mode_III_slip_weakening "
              << "output folder: " << output_folder << "\n"
              << "nb_elements alog x: " << nex << "\n"
              << "nb_time_steps: " << nb_time_steps << "\n"
              << "griffith crack length: " << G_length << "\n"
              << "reference number of elements: " << n_ele_ind << '\n'
              << "dt: " << 0.2 * dx / cs << '\n'
              << std::endl;

    /* -------------------------------------------------------------------------- */

    UInt t = 0;
    UInt x_tip = 0;
    UInt x_lap = 0.05 * nex;

    SpectralModel *model;

    model = new SpectralModel({nex, nez}, 0, {dom_sizex, dom_sizez},
                              nu, E, cs, tcut,
                              sim_name, output_folder);

    Real beta = 0.2;
    // SimulationDriver sim_driver(*model, beta=beta);
    // SimulationDriver sim_driver(*model);
    model->initModel();

    DataRegister::registerParameter("critical_normal_opening", crit_n_open);
    DataRegister::registerParameter("critical_shear_opening", crit_s_open);
    DataRegister::registerParameter("static_friction_coefficient", mus);
    DataRegister::registerParameter("dynamic_friction_coefficient", mud);
    DataRegister::registerParameter("uniform_contact_pressure", load_nor);
    

    Interfacer<_cohesive_coulomb> interfacer(*model);
    interfacer.createUniformInterface();

    interfacer.createThroughCrack((dom_sizex - crack_size) / 2., (dom_sizex + crack_size) / 2.);

    CohesiveLawCoulomb &cohesive_law = dynamic_cast<CohesiveLawCoulomb &>((model->getInterfaceLaw()));

    if(exp) cohesive_law.initExpFormulation();
    else   cohesive_law.initStandardFormulation();

    // cohesive_law.initRegularFormulation();
    //  sim_driver.initConstantLoading(load, psi, phi);
    Real x_mid = nex / 2.0;
    Real z_mid = nez / 2.0;
    Real sint, cost, d;
    // initialize load
    std::vector<Real> initload(nex * nez * 3, 0.0);
    std::vector<Real> load_actu(nex * nez * 3, 0.0);
    for (UInt x = 0; x < nex; x++)
    {
        for (UInt z = 0; z < nez; z++)
        {
            initload[x * 3 + 3 * nex * z + 2] = 0;
            initload[x * 3 + 3 * nex * z + 1] = -load_nor;
            initload[x * 3 + 3 * nex * z + 0] = load_shr;
        }
    }
    // actualload = fluidoverpressure(initload, nex, delta_p_star, alpha, dx, 0);
    model->setLoadingFromVector(initload);
    model->initInterfaceFields();
    /* -------------------------------------------------------------------------- */
    // Set-up simulation  outputs

    // DataDumper dumper(*model);
    DataDumper dumper(*model);

    //dumper.initVectorDumper("ST_Diagram_top_z_velo.cra", _top_velocities, 2, 1.0, 1, 0, _binary);
    //dumper.initVectorDumper("ST_Diagram_top_z_displ.cra", _top_displacements, 2, 1.0, 1, 0, _binary);
    dumper.initVectorDumper("ST_Diagram_top_x_velo.cra", _top_velocities, 1, 1.0, 1, 0, _binary);
    dumper.initVectorDumper("ST_Diagram_top_x_displ.cra", _top_displacements, 1, 1.0, 1, 0, _binary);
    //dumper.initVectorDumper("ST_Diagram_shear_stress.cra", _interface_tractions, 0, 1.0, 1, 0, _binary);
    dumper.initDumper("ST_Diagram_normal_stress.cra", _interface_tractions, 1.0, 1, 0,_binary);
    dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1.0, 1, 0, _binary);


    /* -------------------------------------------------------------------------- */

    // sim_driver.launchCrack(dom_sizex/2.,1.75*G_length,0.075,false);

    // DataRegister::restart_dir = "restart_files/";
    // model->restartModel();
    UInt nb_dumps = 2000;
    UInt nb_t = nb_time_steps / nb_dumps;
    std::ofstream outputFile(output_folder + "ST_cra_tip.cra");
    while ((t < nb_time_steps) && (x_tip < 0.9 * nex))
    {

        // model->pauseModel();
        // break;

        model->updateDisplacements();
        model->fftOnDisplacements();
        model->computeStress();
        model->computeInterfaceFields();
        model->increaseTimeStep();

        x_tip = model->getCohesiveTipPosition(nex / 2, nex, nez/2);

        if (t % nb_t == 0)
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
        // update loading case
        Real time = t * beta * dx / cs;
        load_actu = fluidoverpressure(initload, nex, nez, delta_p_star, alpha, dx, time);

        model->setLoadingFromVector(load_actu);
    }
    outputFile.close();
    model->pauseModel();

    // delete model;
    return 0;
}