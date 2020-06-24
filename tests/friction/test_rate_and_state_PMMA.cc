/**
 * @file   test_rate_and_state_PMMA.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed May  3 17:17:10 2017
 *
 * @brief  Testing onset of sliding along rate and state interface made of PMMA
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
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

  Real k=1.0;
  Real epsilon=1e-7;
  Real D=5e-7;
  Real f_0=0.285;
  Real a=0.005;
  Real b=0.075;
  Real v_star=1e-7;
  Real phi_star=0.00033;
  Real load = 0.36e6;

  // Geometry description 
  Real nu =  0.33;
  Real E = 7.98e8;
  //  Real cs = 1200;
  Real cs = 500;
  // Cut of the loaded material kernels
  UInt tcut = 100;
  
  // Loading case
  Real psi = 90;
  Real phi = 90;
  Real sigma_0=1e6;

  UInt nb_time_steps = 0;
  UInt nb_elements = 512;   
  Real dom_size = 0.5;
    
  /* -------------------------------------------------------------------------- */

  SpectralModel model({nb_elements,1}, nb_time_steps, {dom_size,0.}, nu, nu, 
		      E, E, cs, cs, tcut, tcut,
		      "Testing rate and state friction along PMMA");
  
  model.initModel();

  DataRegister::registerParameter("D_hom",D);
  DataRegister::registerParameter("f_0_hom",f_0);
  DataRegister::registerParameter("a_hom",a);
  DataRegister::registerParameter("b_hom",b);
  DataRegister::registerParameter("v_star_hom",v_star);
  DataRegister::registerParameter("phi_star_hom",phi_star);
  DataRegister::registerParameter("sigma_0",sigma_0);
  
  Interfacer<_rate_and_state> interfacer(model);
  interfacer.createUniformInterface();

  model.setLoadingCase(load, psi, phi);
  model.updateLoads();

  UInt t = 0;

  const CrackProfile * shear_velo_jump = model.readData(_shear_velocity_jumps);
  
  RateAndStateLaw& r_and_s = dynamic_cast<RateAndStateLaw&>((model.getInterfaceLaw()));

  r_and_s.setVelocityPredictor({0.,0.,3e-4});
  
  while (t<20001) {

    model.updateDisplacements();
    model.fftOnDisplacements();
    model.computeStress();

    if (t==0) 
        model.initInterfaceFields();
    else
      model.computeInterfaceFields();

    model.increaseTimeStep();

    if(t==5)
      r_and_s.perturbState(epsilon,k);
    
      if (t%50==0){
      Real v_max = shear_velo_jump->getMaxValue();
      const Real v_av = (*shear_velo_jump)[0];
      std::cout << "Simulation at t " << t << " = " << model.getTime()*dom_size/cs << " [sec]"<< std::endl
		<< "-> v_max-v_av = " << (v_max-v_av)*cs << std::endl;
      }
    ++t;
  }
   
  return 0;
}
