/**
 * @file   travelling_pulse.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Jan 17 11:58:19 2018
 *
 * @brief  Simulation in the framework of travelling pulse solutions
 *
 * @section LICENSE
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
#include "simulation_driver.hh"
#include "interfacer.hh"
#include "cohesive_law.hh"
#include "rate_and_state_law.hh"
#include "regularized_coulomb_law.hh"
#include "data_dumper.hh"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <string>

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]){

  // Note : Construct the pre-integrated material kernels before running this simulation
  // Use "invert_serial.f" to construct kernel files

  std::cout << "./travelling_pulse [ <output_folder_name>='./' <nb_elements>=1024 <dom_size>=20 <patch_size/dom_size>=0.5 <load x tau_m>=1.044 <vm>=0.01 ]"
	    << std::endl;

  std::string sim_name = "searching pulse solutions";
  std::string output_folder = "./";
  Real load = 0.36e6; //0.575e6;
  Real dom_size = 20.;
  UInt nb_elements = 1024;
  Real vm = 0.01;
  Real patch_size = 0.5;
  
  if(argc > 1)
    output_folder = argv[1];
  if(argc > 2)
    nb_elements = std::atof(argv[2]);  
  if(argc > 3)
    dom_size = std::atof(argv[3]);    
  if(argc > 4)
    patch_size = std::atof(argv[4]);    
  if(argc > 5)
    load = std::atof(argv[5])*0.34649e6;      
  if(argc > 6)
    vm = std::atof(argv[6]);    

  // Geometry description 
  Real mu = 9e9;
  Real nu =  0.33;
  Real E = 2*mu*(1+nu);
  //  Real cs = 1200;
  Real cs = 2738.6;

  //Perturbation patch
  std::vector<UInt> limits = {(UInt)((0.5-0.5*patch_size)*nb_elements),(UInt)((0.5+0.5*patch_size)*nb_elements)};
  
  Real std_dev = patch_size*nb_elements/5.;

  // Cut of the loaded material kernels
  UInt tcut = 100;

  UInt nb_simu_steps = 15*nb_elements; //rule of thumb
  
  // Loading case
  Real psi = 90;
  Real phi = 90;
      
  /* -------------------------------------------------------------------------- */

  SpectralModel model(nb_elements, 0, dom_size,
		      nu, E, cs, tcut, sim_name, output_folder);
  
  model.initModel(0.1);
  model.readInputFile("input_pulse_PMMA.dat");

  Real v_predictor = DataRegister::getParameter("v0_predictor");
  UInt t_char = (UInt)(DataRegister::getParameter("dumping_every_t"));

  std::cout << "Input parameters summary : "<< std::endl
	    << "Output directory : " << output_folder << std::endl
	    << "Standard deviation : " << patch_size*dom_size/5. << std::endl
	    << "Patch size : " << patch_size << std::endl
	    << "Number of time steps : " << nb_simu_steps << std::endl
	    << "Patch velocity vm : " << vm << std::endl;

  Interfacer<_regularized_rate_and_state> interfacer(model);
  interfacer.createUniformInterface();

  model.setLoadingCase(load, psi, phi);
  model.updateLoads();
    
  DataDumper dumper(model);
  dumper.initVectorDumper("ST_Diagram_top_z_velo.cra", _top_velocities, 2, 1.0, 1, 0, _binary);
  dumper.initVectorDumper("ST_Diagram_top_z_displ.cra", _top_displacements, 2, 1.0, 1, 0, _binary);
  dumper.initVectorDumper("ST_Diagram_shear_stress.cra", _interface_tractions, 2, 1.0, 1, 0, _binary);
  dumper.initDumper("ST_Diagram_state_variable.cra", _state_variable, 1.0, 1, 0, _binary);
  dumper.initDumper("ST_Diagram_fric_coef.cra", _friction_coefficient, 1.0, 1, 0, _binary);
  
  UInt t = 0;

  const CrackProfile * shear_velo_jump = model.readData(_shear_velocity_jumps);

  RateAndStateLaw& r_and_s = dynamic_cast<RateAndStateLaw&>((model.getInterfaceLaw()));
  
  r_and_s.setVelocityPredictor({0.,0.,v_predictor});
  
  Real v_av,v_5;
  bool dynamic = false;

  //  while (v_25*cs<v_min_fric) {
  while (t<nb_simu_steps) {

    model.updateDisplacements();
    model.fftOnDisplacements();
    model.computeStress();

    if (t==0)
      model.initInterfaceFields();
    else
      model.computeInterfaceFields();

    model.increaseTimeStep();

    if (t==5)
      r_and_s.insertGaussianPerturbation(std_dev,vm);
      //r_and_s->insertPerturbationPatch(limits,vm);

    if ((t%t_char==0)||((dynamic)&&(t%5==0))||(t==5))
      dumper.dumpAll();
    
    if (t%(3*t_char)==0){

      v_av = (*shear_velo_jump)[0];
      v_5 = (*shear_velo_jump)[(int)(0.5*nb_elements)];

      std::cout << "Simulation at t " << t << " = " << model.getTime()*dom_size/cs << " [sec]"<< std::endl
		<< "sliding velocity at the center point: " << v_5*cs << std::endl
		<< "sliding velocity at an edge point: " << v_av*cs << std::endl;
      }
    ++t;
  }
  return 0;
}
 
