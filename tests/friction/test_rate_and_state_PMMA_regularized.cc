/**
 * @file   test_rate_and_state_PMMA_regularized.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Tue Jan 16 18:52:03 2018
 *
 * @brief  Testing the regularized law
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
  Real load = 1.005*0.34475668e6;
  Real theta = 5./18.;
  Real xi = 0.0055;
  Real v0 = 1e-9;
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
  UInt nb_elements = 2048;   
  Real dom_size = 2.0;
    
  /* -------------------------------------------------------------------------- */

  SpectralModel model(nb_elements, nb_time_steps, dom_size,
		      nu, E, cs, tcut, 
		      "Testing rate and state friction along PMMA");
  
  model.initModel();

  DataRegister::registerParameter("D_hom",D);
  DataRegister::registerParameter("f_0_hom",f_0);
  DataRegister::registerParameter("a_hom",a);
  DataRegister::registerParameter("b_hom",b);
  DataRegister::registerParameter("v_star_hom",v_star);
  DataRegister::registerParameter("phi_star_hom",phi_star);
  DataRegister::registerParameter("sigma_0",sigma_0);
  DataRegister::registerParameter("theta",theta);
  DataRegister::registerParameter("xi",xi);
  DataRegister::registerParameter("v0",v0);
  
  Interfacer<_regularized_rate_and_state> interfacer(model);
  RateAndStateLaw& r_and_s = dynamic_cast<RateAndStateLaw&>((model.getInterfaceLaw()));
  r_and_s.initRegularizedStateEvolution(DataRegister::getParameter<Real>("v0"));
  interfacer.createUniformInterface();

  model.setLoadingCase(load, psi, phi);
  model.updateLoads();

  UInt t = 0;

  const CrackProfile * shear_velo_jump = model.readData(_shear_velocity_jumps);
  const std::vector<Real> * state = model.readData(_state_variable);
    
  
  r_and_s.setVelocityPredictor({0.,0.,3e-4});
  
  Real v_max=0.;
  Real v_av,v_25,state_75,state_ss; 
  
  while (t<5001) {

    model.updateDisplacements();
    model.fftOnDisplacements();
    model.computeStress();

    if (t==0) {
      model.initInterfaceFields();
      state_ss = (*state)[(int)(0.75*nb_elements)];
    }
    else
      model.computeInterfaceFields();
    
    model.increaseTimeStep();

    if(t==5)
      r_and_s.perturbState(epsilon,k);
    
      if (t%50==0){
	v_max = shear_velo_jump->getMaxValue();
	v_av = (*shear_velo_jump)[0];
	v_25 = (*shear_velo_jump)[(int)(0.25*nb_elements)];
	state_75 = (*state)[(int)(0.75*nb_elements)];
	
	std::cout << "Simulation at t " << t << " = " << model.getTime() << " [sec]"<< std::endl
		  << "sliding velocity at the observation point: " << v_25*cs << std::endl
		  << "state variable at the expected nucleation point: " << state_75/state_ss << std::endl
		  << "-> v_max-v_av = " << (v_max-v_av)*cs << std::endl;
      }
      ++t;
  }
   
  return 0;
}
