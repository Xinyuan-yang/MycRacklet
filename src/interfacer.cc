/* -------------------------------------------------------------------------- */
#include "interfacer.hh"
/* -------------------------------------------------------------------------- */
template void Interfacer<_rate_and_state>::createHomogeneousRateandStateIntfc();
template void Interfacer<_weakening_rate_and_state>::createHomogeneousRateandStateIntfc();
template void Interfacer<_regularized_rate_and_state>::createHomogeneousRateandStateIntfc();
template void Interfacer<_regularized_weakening_rate_and_state>::createHomogeneousRateandStateIntfc();

template<FractureLawType F>
void Interfacer<F>::createHomogeneousRateandStateIntfc() {
  
  Real D_value = DataRegister::getParameter("D_hom");
  Real f_0_value = DataRegister::getParameter("f_0_hom");
  Real a_value = DataRegister::getParameter("a_hom");
  Real b_value = DataRegister::getParameter("b_hom");
  Real v_star_value = DataRegister::getParameter("v_star_hom");
  Real phi_star_value = DataRegister::getParameter("phi_star_hom");
  
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

  Real sigma_0 = DataRegister::getParameter("sigma_0");
  
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
