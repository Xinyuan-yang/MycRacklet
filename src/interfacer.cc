/**
 * @file   interfacer.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Sep 11 13:17:04 2014
 *
 * @brief Implementation of the Interfacer class
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
 */
/* -------------------------------------------------------------------------- */
#include "interfacer.hh"
#include <iostream>
#include <fstream>
/* -------------------------------------------------------------------------- */

// Homogeneous rate and state friction interfaces

template void Interfacer<_rate_and_state>::createHomogeneousRateandStateIntfc();
template void Interfacer<_weakening_rate_and_state>::createHomogeneousRateandStateIntfc();
template void Interfacer<_regularized_rate_and_state>::createHomogeneousRateandStateIntfc();
template void Interfacer<_regularized_weakening_rate_and_state>::createHomogeneousRateandStateIntfc();

// Heterogeneous rate and state friction interfaces

template void Interfacer<_rate_and_state>::createHeterogeneousRateandStateIntfc(std::vector<Real> vec_D, std::vector<Real> vec_f0, std::vector<Real> vec_a, std::vector<Real> vec_b, std::vector<Real> vec_v_star, std::vector<Real> vec_phi_star);
template void Interfacer<_weakening_rate_and_state>::createHeterogeneousRateandStateIntfc(std::vector<Real> vec_D, std::vector<Real> vec_f0, std::vector<Real> vec_a, std::vector<Real> vec_b, std::vector<Real> vec_v_star, std::vector<Real> vec_phi_star);
template void Interfacer<_regularized_rate_and_state>::createHeterogeneousRateandStateIntfc(std::vector<Real> vec_D, std::vector<Real> vec_f0, std::vector<Real> vec_a, std::vector<Real> vec_b, std::vector<Real> vec_v_star, std::vector<Real> vec_phi_star);
template void Interfacer<_regularized_weakening_rate_and_state>::createHeterogeneousRateandStateIntfc(std::vector<Real> vec_D, std::vector<Real> vec_f0, std::vector<Real> vec_a, std::vector<Real> vec_b, std::vector<Real> vec_v_star, std::vector<Real> vec_phi_star);

// Homegeneous cohesive coulomb interfaces
template void Interfacer<_cohesive_coulomb>::createHomogeneousCoulombIntfc();
template void Interfacer<_dual_cohesive_coulomb>::createHomogeneousCoulombIntfc();

template<FractureLawType F>
void Interfacer<F>::createHomogeneousCoulombIntfc() {
  Real crit_nor_opening = DataRegister::getParameter<Real>("critical_normal_opening");
  Real crit_shr_opening = DataRegister::getParameter<Real>("critical_shear_opening");
  Real coef_static = DataRegister::getParameter<Real>("static_friction_coefficient");
  Real coef_dynamic = DataRegister::getParameter<Real>("dynamic_friction_coefficient");
  Real sigma_0 = DataRegister::getParameter<Real>("uniform_contact_pressure");

  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * cf_s = datas[_static_friction_coefficient];
  std::vector<Real> * cf_d = datas[_dynamic_friction_coefficient];
  std::vector<Real> * fric_strength = datas[_frictional_strength];

  std::fill(crit_n_open->begin(), crit_n_open->end(), crit_nor_opening);
  std::fill(crit_s_open->begin(), crit_s_open->end(), crit_shr_opening);
  std::fill(cf_s->begin(), cf_s->end(), coef_static);
  std::fill(cf_d->begin(), cf_d->end(), coef_dynamic);
  std::fill(fric_strength->begin(), fric_strength->end(), coef_static*sigma_0);

  // For DualCohesiveCoulombFormulation
  if (DataRegister::hasParameter("critical_int_opening")) {
    Real crit_int_opening = DataRegister::getParameter<Real>("critical_int_opening");  
    std::vector<Real> * crit_i_open = datas[_critical_int_opening];
    std::fill(crit_i_open->begin(), crit_i_open->end(), crit_int_opening);
  }
  if (DataRegister::hasParameter("int_friction_coefficient")) {
    Real coeff_int = DataRegister::getParameter<Real>("int_friction_coefficient");  
    std::vector<Real> * cf_int = datas[_critical_int_opening];
    std::fill(cf_int->begin(), cf_int->end(), coeff_int);
  }

  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " UNIFORM INTERFACE " << std::endl
        << "* Critical normal opening: " << crit_nor_opening << std::endl
	      << "* Critical shear opening: " << crit_shr_opening << std::endl
	      << "* Static friction coefficient: " << coef_static << std::endl
	      << "* Dynamic friction coefficient: " << coef_dynamic << std::endl
	      << "* Uniform contact pressure: " << sigma_0 << std::endl
	      << std::endl;	

  out_parameters << "delta_c_nor " << crit_nor_opening << std::endl
		    << "delta_c_shr " << crit_shr_opening << std::endl
        << "cf_s " << coef_static << std::endl
	      << "cf_d " << coef_dynamic << std::endl
	      << "sigma_0 " << sigma_0 << std::endl;  
}

template<FractureLawType F>
void Interfacer<F>::createHomogeneousRateandStateIntfc() {
  
  Real D_value = DataRegister::getParameter<Real>("D_hom");
  Real f_0_value = DataRegister::getParameter<Real>("f_0_hom");
  Real a_value = DataRegister::getParameter<Real>("a_hom");
  Real b_value = DataRegister::getParameter<Real>("b_hom");
  Real v_star_value = DataRegister::getParameter<Real>("v_star_hom");
  Real phi_star_value = DataRegister::getParameter<Real>("phi_star_hom");
  
  std::vector<Real> * D = datas[_rands_D];
  std::vector<Real> * f_0 = datas[_rands_f_0];
  std::vector<Real> * a = datas[_rands_a];
  std::vector<Real> * b = datas[_rands_b];
  std::vector<Real> * v_star = datas[_rands_v_star];
  std::vector<Real> * phi_star = datas[_rands_phi_star];

  std::fill(D->begin(),D->end(), D_value);
  std::fill(f_0->begin(),f_0->end(), f_0_value);
  std::fill(a->begin(),a->end(), a_value);
  std::fill(b->begin(),b->end(), b_value);
  std::fill(v_star->begin(),v_star->end(), v_star_value);
  std::fill(phi_star->begin(),phi_star->end(), phi_star_value);

  Real sigma_0 = DataRegister::getParameter<Real>("sigma_0");
  
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " HOMOGENEOUS RATE AND STATE FRICTIONAL INTERFACE: " << std::endl
	      << "* D: " << D_value << std::endl    
	      << "* f_0: " << f_0_value << std::endl    
	      << "* a: " << a_value << std::endl    
	      << "* b: " << b_value << std::endl  
	      << "* v_star: " << v_star_value << std::endl    
	      << "* phi_star: " << phi_star_value << std::endl    
	      << "* Uniform contact pressure: " << sigma_0 << std::endl
	      << std::endl;	
  out_parameters << "D " << D_value << std::endl
		 << "f_0 " << f_0_value << std::endl
		 << "a " << a_value << std::endl
		 << "b " << b_value << std::endl
		 << "v_star " << v_star_value << std::endl
		 << "phi_star " << phi_star_value << std::endl
		 << "sigma_0 " << sigma_0 << std::endl;
}


template<FractureLawType F>
void Interfacer<F>::createHeterogeneousRateandStateIntfc(std::vector<Real> vec_D, std::vector<Real> vec_f0, std::vector<Real> vec_a, std::vector<Real> vec_b, std::vector<Real> vec_v_star, std::vector<Real> vec_phi_star) {
    
  std::vector<Real> * D = datas[_rands_D];
  std::vector<Real> * f_0 = datas[_rands_f_0];
  std::vector<Real> * a = datas[_rands_a];
  std::vector<Real> * b = datas[_rands_b];
  std::vector<Real> * v_star = datas[_rands_v_star];
  std::vector<Real> * phi_star = datas[_rands_phi_star];

  for (UInt x = 0; x < n_ele[0]; ++x) {
    for (UInt z = 0; z < n_ele[1]; ++z ) {
      (*D)[x+z*n_ele[0]] = vec_D[x+z*n_ele[0]];
      (*f_0)[x+z*n_ele[0]] = vec_f0[x+z*n_ele[0]];
      (*a)[x+z*n_ele[0]] = vec_a[x+z*n_ele[0]];
      (*b)[x+z*n_ele[0]] = vec_b[x+z*n_ele[0]];
      (*v_star)[x+z*n_ele[0]] = vec_v_star[x+z*n_ele[0]];
      (*phi_star)[x+z*n_ele[0]] = vec_phi_star[x+z*n_ele[0]];
    }
  }
  
  Real sigma_0 = DataRegister::getParameter<Real>("sigma_0");
  
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " HETEROGENEOUS RATE AND STATE FRICTIONAL INTERFACE!! " << std::endl
	      << "* Uniform contact pressure: " << sigma_0 << std::endl
	      << std::endl;	
}
