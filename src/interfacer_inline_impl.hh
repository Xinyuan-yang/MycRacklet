/**
 * @file   interfacer_inline_impl.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Sep 11 13:17:04 2014
 *
 * @brief  Implementation of the inline functions of the Interfacer class
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
#include "cohesive_law.hh"
#include "cohesive_law_viscoelastic.hh"
#include "cohesive_law_all.hh"
#include "cohesive_law_coulomb.hh"
#include <random>
#include <chrono>
/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::initInterfaceLaw() {
  interface_law = std::make_shared<CohesiveLaw>(); 
}

template<>
inline void Interfacer<_coupled_cohesive>::initInterfaceLaw() {
  interface_law = std::make_shared<CohesiveLawAll>(); 
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_viscoelastic_coupled_cohesive>::initInterfaceLaw() {
  interface_law = std::make_shared<CohesiveLawViscoelastic>(); 
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_cohesive_coulomb>::initInterfaceLaw() {
  interface_law = std::make_shared<CohesiveLawCoulomb>(); 
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_dual_cohesive_coulomb>::initInterfaceLaw() {
  interface_law = std::make_shared<CohesiveLawCoulomb>(); 
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_rate_and_state>::initInterfaceLaw() {
  interface_law = std::make_shared<RateAndStateLaw>();
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_weakening_rate_and_state>::initInterfaceLaw() {
  interface_law = std::make_shared<RateAndStateLaw>();
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_regularized_rate_and_state>::initInterfaceLaw() {
  interface_law = std::make_shared<RateAndStateLaw>();
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_regularized_weakening_rate_and_state>::initInterfaceLaw() {
  interface_law = std::make_shared<RateAndStateLaw>();
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_rate_and_state>::createUniformInterface() {

  createHomogeneousRateandStateIntfc();
  std::shared_ptr<RateAndStateLaw> r_and_s = std::dynamic_pointer_cast<RateAndStateLaw>(interface_law);
  r_and_s->initStandardFormulation();
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_weakening_rate_and_state>::createUniformInterface() {

  createHomogeneousRateandStateIntfc();
  std::shared_ptr<RateAndStateLaw> r_and_s = std::dynamic_pointer_cast<RateAndStateLaw>(interface_law);
  r_and_s->initVelocityWeakeningFormulation();
  
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_regularized_rate_and_state>::createUniformInterface() {

  Real v0 = DataRegister::getParameter<Real>("v0");
   
  createHomogeneousRateandStateIntfc();
  std::shared_ptr<RateAndStateLaw> r_and_s = std::dynamic_pointer_cast<RateAndStateLaw>(interface_law);
  r_and_s->initRegularizedFormulation(v0);
  
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_regularized_rate_and_state>::createHeterogeneousInterface(std::vector<Real> vec_D, std::vector<Real> vec_f0, std::vector<Real> vec_a, std::vector<Real> vec_b, std::vector<Real> vec_v_star, std::vector<Real> vec_phi_star) {

  Real v0 = DataRegister::getParameter<Real>("v0");
   
  createHeterogeneousRateandStateIntfc(vec_D,vec_f0,vec_a,vec_b,vec_v_star,vec_phi_star);
  std::shared_ptr<RateAndStateLaw> r_and_s = std::dynamic_pointer_cast<RateAndStateLaw>(interface_law);
  r_and_s->initRegularizedFormulation(v0);
  
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_regularized_weakening_rate_and_state>::createUniformInterface() {

  Real v0 = DataRegister::getParameter<Real>("v0");
   
  createHomogeneousRateandStateIntfc();
  std::shared_ptr<RateAndStateLaw> r_and_s = std::dynamic_pointer_cast<RateAndStateLaw>(interface_law);
  r_and_s->initRegularizedWeakeningFormulation(v0);
  
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_regularized_weakening_rate_and_state>::createHeterogeneousInterface(std::vector<Real> vec_D, std::vector<Real> vec_f0, std::vector<Real> vec_a, std::vector<Real> vec_b, std::vector<Real> vec_v_star, std::vector<Real> vec_phi_star) {

  Real v0 = DataRegister::getParameter<Real>("v0");
   
  createHeterogeneousRateandStateIntfc(vec_D,vec_f0,vec_a,vec_b,vec_v_star,vec_phi_star);
  std::shared_ptr<RateAndStateLaw> r_and_s = std::dynamic_pointer_cast<RateAndStateLaw>(interface_law);
  r_and_s->initRegularizedWeakeningFormulation(v0);
  
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_cohesive_coulomb>::createUniformInterface() {

  createHomogeneousCoulombIntfc();
  std::shared_ptr<CohesiveLawCoulomb> cohesive_law = std::dynamic_pointer_cast<CohesiveLawCoulomb>(interface_law);
  cohesive_law->initStandardFormulation();
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_dual_cohesive_coulomb>::createUniformInterface() {

  createHomogeneousCoulombIntfc();
  std::shared_ptr<CohesiveLawCoulomb> cohesive_law = std::dynamic_pointer_cast<CohesiveLawCoulomb>(interface_law);
  cohesive_law->initDualFormulation();
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createUniformInterface() {

  Real crit_nor_opening = DataRegister::getParameter<Real>("critical_normal_opening");
  Real crit_shr_opening = DataRegister::getParameter<Real>("critical_shear_opening");
  Real max_nor_strength = DataRegister::getParameter<Real>("max_normal_strength");
  Real max_shr_strength = DataRegister::getParameter<Real>("max_shear_strength");
  
  Real res_nor_strength = 0;
  Real res_shr_strength = 0;
  if(DataRegister::hasParameter("res_normal_strength")){
    res_nor_strength = DataRegister::getParameter<Real>("res_normal_strength");
  }
  if(DataRegister::hasParameter("res_shear_strength")){
    res_shr_strength = DataRegister::getParameter<Real>("res_shear_strength");
  }
  
  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];
  
  std::fill(nor_strength->begin(), nor_strength->end(), max_nor_strength);
  std::fill(shr_strength->begin(), shr_strength->end(), max_shr_strength);
  std::fill(crit_n_open->begin(), crit_n_open->end(), crit_nor_opening);
  std::fill(crit_s_open->begin(), crit_s_open->end(), crit_shr_opening);
  std::fill(res_n_strength->begin(), res_n_strength->end(), res_nor_strength);
  std::fill(res_s_strength->begin(), res_s_strength->end(), res_shr_strength);
  
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " UNIFORM INTERFACE " << std::endl
	      << "* Critical normal opening: " << crit_nor_opening << std::endl
	      << "* Critical shear opening: " << crit_shr_opening << std::endl
	      << "* Maximal normal strength: " << max_nor_strength << std::endl
	      << "* Maximal shear strength: " << max_shr_strength << std::endl
	      << "* Residual normal strength: " << res_nor_strength << std::endl
	      << "* Residual shear strength: " << res_shr_strength << std::endl
	      << std::endl;	

  out_parameters << "delta_c_nor " << crit_nor_opening << std::endl
		 << "delta_c_shr " << crit_shr_opening << std::endl
		 << "tau_max_nor " << max_nor_strength << std::endl
		 << "tau_max_shr " << max_shr_strength << std::endl
    		 << "tau_res_nor " << res_nor_strength << std::endl
		 << "tau_res_shr " << res_shr_strength << std::endl;  
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_coupled_cohesive>::createUniformInterface() {

  Real crit_nor_opening = DataRegister::getParameter<Real>("critical_normal_opening");
  Real crit_shr_opening = DataRegister::getParameter<Real>("critical_shear_opening");
  Real max_nor_strength = DataRegister::getParameter<Real>("max_normal_strength");
  Real max_shr_strength = DataRegister::getParameter<Real>("max_shear_strength");
  
  Real res_nor_strength = 0;
  Real res_shr_strength = 0;
  if(DataRegister::hasParameter("res_normal_strength")){
    res_nor_strength = DataRegister::getParameter<Real>("res_normal_strength");
  }
  if(DataRegister::hasParameter("res_shear_strength")){
    res_shr_strength = DataRegister::getParameter<Real>("res_shear_strength");
  }
  
  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];
  
  std::fill(nor_strength->begin(), nor_strength->end(), max_nor_strength);
  std::fill(shr_strength->begin(), shr_strength->end(), max_shr_strength);
  std::fill(crit_n_open->begin(), crit_n_open->end(), crit_nor_opening);
  std::fill(crit_s_open->begin(), crit_s_open->end(), crit_shr_opening);
  std::fill(res_n_strength->begin(), res_n_strength->end(), res_nor_strength);
  std::fill(res_s_strength->begin(), res_s_strength->end(), res_shr_strength);
  
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " UNIFORM INTERFACE " << std::endl
	      << "* Critical normal opening: " << crit_nor_opening << std::endl
	      << "* Critical shear opening: " << crit_shr_opening << std::endl
	      << "* Maximal normal strength: " << max_nor_strength << std::endl
	      << "* Maximal shear strength: " << max_shr_strength << std::endl
	      << "* Residual normal strength: " << res_nor_strength << std::endl
	      << "* Residual shear strength: " << res_shr_strength << std::endl
	      << std::endl;	

  out_parameters << "delta_c_nor " << crit_nor_opening << std::endl
		 << "delta_c_shr " << crit_shr_opening << std::endl
		 << "tau_max_nor " << max_nor_strength << std::endl
		 << "tau_max_shr " << max_shr_strength << std::endl
    		 << "tau_res_nor " << res_nor_strength << std::endl
		 << "tau_res_shr " << res_shr_strength << std::endl;  
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::insertPatternfromFile(std::string filename, UInt origin) {

  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];

  std::fill(res_n_strength->begin(),  res_n_strength->end(), 0);
  std::fill(res_s_strength->begin(),  res_s_strength->end(), 0);
  
  std::ifstream file;
  std::string line;
  file.open(filename);
  Real ratio;
  
  for (UInt x = origin; x < n_ele[0]; ++x) {
    if(file.eof())
      break;
    std::getline(file,line);
    std::stringstream sstr(line);
    for (UInt z = 0; z < n_ele[1]; ++z ) {
      sstr >> ratio;
      (*shr_strength)[x+z*n_ele[0]]*= ratio;
      (*nor_strength)[x+z*n_ele[0]]*= ratio;
    }
  }
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " TOUGHNESS PATTERN INSERTED FROM FILE " << std::endl
	      << "* Filename: " << filename << std::endl
	      << std::endl;	
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::insertPatternfromFile(std::string file_strength, std::string file_opening, std::string file_residual) {

  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];

  std::ifstream file_str;
  std::string line_str;
  file_str.open(file_strength);
  Real strength;

  std::ifstream file_op;
  std::string line_op;
  file_op.open(file_opening);
  Real opening;
  
  for (UInt x = 0; x < n_ele[0]; ++x) {
    if(file_str.eof()||file_op.eof())
      break;

    // Read from the strength file
    std::getline(file_str,line_str);
    std::stringstream sstr_str(line_str);

    // Read from the opening file
    std::getline(file_op,line_op);
    std::stringstream sstr_op(line_op);

    for (UInt z = 0; z < n_ele[1]; ++z ) {

      sstr_str >> strength;
      sstr_op >> opening;
      
      (*shr_strength)[x+z*n_ele[0]] = strength;
      (*nor_strength)[x+z*n_ele[0]] = strength;
      (*crit_s_open)[x+z*n_ele[0]] = opening;
      (*crit_n_open)[x+z*n_ele[0]] = opening;
      
      if(strength==0) {
	(*ind_crack)[x+z*n_ele[0]] = 0;
      }
    }
  }
  
  if (file_residual == "None"){
    std::fill(res_n_strength->begin(),  res_n_strength->end(), 0);
    std::fill(res_s_strength->begin(),  res_s_strength->end(), 0);
  }
  else {
    std::ifstream file_res;
    std::string line_res;
    file_res.open(file_residual);
    Real residual;

    for (UInt x = 0; x < n_ele[0]; ++x) {
      if(file_res.eof())
	break;

      // Read from the residual file
      std::getline(file_res,line_res);
      std::stringstream sstr_res(line_res);

      for (UInt z = 0; z < n_ele[1]; ++z ) {
	
	sstr_res >> residual;
	(*res_n_strength)[x+z*n_ele[0]] = residual;
	(*res_s_strength)[x+z*n_ele[0]] = residual;
	
      }
    }
  }
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " STRENGTH PATTERN INSERTED FROM FILE " << std::endl
	      << "* Filename: " << file_strength << std::endl
	      << " OPENING PATTERN INSERTED FROM FILE " << std::endl
	      << "* Filename: " << file_opening << std::endl
	      << std::endl;	
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createHeterogeneousInterface(std::vector<Real> crit_nor_opening, std::vector<Real> max_nor_strength,std::vector<Real> res_nor_strength, std::vector<Real> crit_shr_opening, std::vector<Real> max_shr_strength, std::vector<Real> res_shr_strength) {

  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];
  
  // Check the length of the shear arguments...

  if ( crit_shr_opening.size() == 0 ){
    crit_shr_opening = crit_nor_opening;
    max_shr_strength = max_nor_strength;
  }
  
  for (UInt x = 0; x < n_ele[0]; ++x) {
    for (UInt z = 0; z < n_ele[1]; ++z ) {

      (*shr_strength)[x+z*n_ele[0]] = max_shr_strength[x+z*n_ele[0]];
      (*nor_strength)[x+z*n_ele[0]] = max_nor_strength[x+z*n_ele[0]];
      (*crit_s_open)[x+z*n_ele[0]] = crit_shr_opening[x+z*n_ele[0]];
      (*crit_n_open)[x+z*n_ele[0]] = crit_nor_opening[x+z*n_ele[0]];
      
      if(max_shr_strength[x+z*n_ele[0]]==0 || max_nor_strength[x+z*n_ele[0]]) {
	(*ind_crack)[x+z*n_ele[0]] = 0;
      }
    }
  }

  // Check the length of the residual normal
  
  if ( res_nor_strength.size() == 0 ){
    std::fill(res_n_strength->begin(),  res_n_strength->end(), 0);
  }
  else{
    for (UInt x = 0; x < n_ele[0]; ++x) {
      for (UInt z = 0; z < n_ele[1]; ++z ) {
	(*res_n_strength)[x+z*n_ele[0]] = res_nor_strength[x+z*n_ele[0]];
      }
    }
  }
  
  // Check the length of the residual shear
  
  if ( res_shr_strength.size() == 0 ){
    std::fill(res_s_strength->begin(),  res_s_strength->end(), 0);
  }
  else{
    for (UInt x = 0; x < n_ele[0]; ++x) {
      for (UInt z = 0; z < n_ele[1]; ++z ) {
	(*res_s_strength)[x+z*n_ele[0]] = res_shr_strength[x+z*n_ele[0]];
      }
    }
  }    
  
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << "Maximum strength and critical opening patterns inserted from vectors."
	      << std::endl
	      << std::endl;
}

/* -------------------------------------------------------------------------- */
template<>inline
void Interfacer<_linear_coupled_cohesive>::createNormalDistributedInterface(Real crit_nor_opening, 
									    Real crit_shr_opening, 
									    Real max_nor_strength, 
									    Real max_shr_strength,
									    Real stddev, Real seed) {
 
  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];

  std::fill(res_n_strength->begin(),  res_n_strength->end(), 0);
  std::fill(res_s_strength->begin(),  res_s_strength->end(), 0);
  
  // construct a trivial random generator engine from a time-based seed:
  //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::normal_distribution<Real> distribution (max_shr_strength,stddev);

  for (UInt i = 0; i < total_n_ele; ++i) {

    Real random_strength = distribution(generator);

    (*nor_strength)[i] = random_strength;
    (*shr_strength)[i] = random_strength;
    (*crit_n_open)[i] = crit_nor_opening;
    (*crit_s_open)[i] = crit_shr_opening;
    (*ind_crack)[i] = 0;
  }
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " HETEROGENEOUS INTERFACE FOLLOWING A NORMAL DISTRIBUTION" << std::endl
	      << "* Standard deviation: " << stddev << std::endl
	      << "* Seed: " << seed << std::endl
    	      << "* Critical normal opening: " << crit_nor_opening << std::endl
	      << "* Critical shear opening: " << crit_shr_opening << std::endl
	      << "* Maximal normal strength: " << max_nor_strength << std::endl
	      << "* Maximal shear strength: " << max_shr_strength << std::endl
	      << std::endl;	
  
  out_parameters << "delta_c_nor " << crit_nor_opening << std::endl
		 << "delta_c_shr " << crit_shr_opening << std::endl
		 << "tau_max_nor " << max_nor_strength << std::endl
		 << "tau_max_shr " << max_shr_strength << std::endl;  
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughArea(Real area_start, Real area_end,
								    UInt cracking_index,
								    Real ratio_max_nor_strength, 
								    Real ratio_max_shr_strength,
								    Real ratio_crit_nor_opening, 
								    Real ratio_crit_shr_opening,
								    bool vrtr) {

  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  
  UInt i_start = (UInt)(area_start/dx[0]);
  UInt i_end = (UInt)(area_end/dx[0]);

  for (UInt ix = i_start; ix < i_end; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt index = ix+iz*n_ele[0];

      (*nor_strength)[index] *= (vrtr+ratio_max_nor_strength);      
      (*shr_strength)[index] *= (vrtr+ratio_max_shr_strength);
      (*crit_n_open)[index] *= (vrtr+ratio_crit_nor_opening);
      (*crit_s_open)[index] *= (vrtr+ratio_crit_shr_opening);
      (*ind_crack)[index] = cracking_index;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_cohesive_coulomb>::createThroughCrack(Real crack_start, Real crack_end) {
  std::vector<Real> * cf_d = datas[_dynamic_friction_coefficient];  
  std::vector<Real> * fric_strength = datas[_frictional_strength]; 
  std::vector<UInt> * ind_crack = datas[_id_crack]; 

  Real sigma_0 = DataRegister::getParameter<Real>("uniform_contact_pressure");

  UInt i_start = (UInt)(crack_start/dx[0]);
  UInt i_end = (UInt)(crack_end/dx[0]);

  for (UInt ix = i_start; ix < i_end; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt index = ix+iz*n_ele[0];

      (*fric_strength)[index] = (*cf_d)[index]*sigma_0;      
      (*ind_crack)[index] = 2;
    }
  }

  Real initial_crack_size = crack_end-crack_start;
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " THROUGH CRACK " << std::endl
	      << "* Initial crack size: " << initial_crack_size << std::endl
	      << "* Crack starts at: " << crack_start << std::endl
	      << "* Crack ends at: " << crack_end << std::endl
	      << "* The Strength has been set to the residual value" << std::endl
	      << std::endl;	

  out_parameters << "a0 " << initial_crack_size << std::endl
		 << "a_l " << crack_start << std::endl
		 << "a_r " << crack_end << std::endl;  
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_dual_cohesive_coulomb>::createThroughCrack(Real crack_start, Real crack_end) {
  std::vector<Real> * cf_d = datas[_dynamic_friction_coefficient];  
  std::vector<Real> * fric_strength = datas[_frictional_strength]; 
  std::vector<UInt> * ind_crack = datas[_id_crack]; 

  Real sigma_0 = DataRegister::getParameter<Real>("uniform_contact_pressure");

  UInt i_start = (UInt)(crack_start/dx[0]);
  UInt i_end = (UInt)(crack_end/dx[0]);

  for (UInt ix = i_start; ix < i_end; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt index = ix+iz*n_ele[0];

      (*fric_strength)[index] = (*cf_d)[index]*sigma_0;    
      (*ind_crack)[index] = 2;
    }
  }

  Real initial_crack_size = crack_end-crack_start;
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " THROUGH CRACK " << std::endl
	      << "* Initial crack size: " << initial_crack_size << std::endl
	      << "* Crack starts at: " << crack_start << std::endl
	      << "* Crack ends at: " << crack_end << std::endl
	      << "* The Strength has been set to the residual value" << std::endl
	      << std::endl;	

  out_parameters << "a0 " << initial_crack_size << std::endl
		 << "a_l " << crack_start << std::endl
		 << "a_r " << crack_end << std::endl;  
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughCrack(Real crack_start, Real crack_end) {
  
  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];

  UInt i_start = (UInt)(crack_start/dx[0]);
  UInt i_end = (UInt)(crack_end/dx[0]);

  for (UInt ix = i_start; ix < i_end; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt index = ix+iz*n_ele[0];

      (*nor_strength)[index] = (*res_n_strength)[index];      
      (*shr_strength)[index] = (*res_s_strength)[index];
      (*ind_crack)[index] = 2;
    }
  }

  Real initial_crack_size = crack_end-crack_start;
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " THROUGH CRACK " << std::endl
	      << "* Initial crack size: " << initial_crack_size << std::endl
	      << "* Crack starts at: " << crack_start << std::endl
	      << "* Crack ends at: " << crack_end << std::endl
	      << "* The Strength has been set to the residual value" << std::endl
	      << std::endl;	

  out_parameters << "a0 " << initial_crack_size << std::endl
		 << "a_l " << crack_start << std::endl
		 << "a_r " << crack_end << std::endl;
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_coupled_cohesive>::createThroughCrack(Real crack_start, Real crack_end) {
  
  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];

  UInt i_start = (UInt)(crack_start/dx[0]);
  UInt i_end = (UInt)(crack_end/dx[0]);

  for (UInt ix = i_start; ix < i_end; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt index = ix+iz*n_ele[0];

      (*nor_strength)[index] = (*res_n_strength)[index];      
      (*shr_strength)[index] = (*res_s_strength)[index];
      (*ind_crack)[index] = 2;
    }
  }

  Real initial_crack_size = crack_end-crack_start;
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " THROUGH CRACK " << std::endl
	      << "* Initial crack size: " << initial_crack_size << std::endl
	      << "* Crack starts at: " << crack_start << std::endl
	      << "* Crack ends at: " << crack_end << std::endl
	      << "* The Strength has been set to the residual value" << std::endl
	      << std::endl;	

  out_parameters << "a0 " << initial_crack_size << std::endl
		 << "a_l " << crack_start << std::endl
		 << "a_r " << crack_end << std::endl;
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughWall(Real wall_start, Real wall_end) {
  
  createThroughArea(wall_start, wall_end, 6, 1000., 1000.);

  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " THROUGH WALL " << std::endl
	      << "* Wall starts at: " << wall_start << std::endl
	      << "* Wall ends at: " << wall_end << std::endl
	      << std::endl;	
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughPolarAsperity(Real position, Real width,
									     Real delta_max_nor_strength, 
									     Real delta_max_shr_strength,
									     Real delta_crit_nor_opening, 
									     Real delta_crit_shr_opening, 
									     bool polarity) {

  if(polarity) {
    createThroughArea(position-width/2,position,5,-delta_max_nor_strength,-delta_max_shr_strength,
		      -delta_crit_nor_opening, -delta_crit_shr_opening, 1);
    createThroughArea(position,position+width/2,4,delta_max_nor_strength,delta_max_shr_strength,
		      delta_crit_nor_opening, delta_crit_shr_opening, 1);
  }
  else {
    createThroughArea(position-width/2,position,4,delta_max_nor_strength,delta_max_shr_strength,
		      delta_crit_nor_opening, delta_crit_shr_opening, 1);
    createThroughArea(position,position+width/2,5,-delta_max_nor_strength,-delta_max_shr_strength,
		      -delta_crit_nor_opening, -delta_crit_shr_opening, 1);
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline UInt Interfacer<_linear_coupled_cohesive>::createThroughMultiPolAsperity(Real start, Real end,
										Real number,
										Real delta_max_nor_strength, 
										Real delta_max_shr_strength,
										Real delta_crit_nor_opening, 
										Real delta_crit_shr_opening, 
										bool polarity) {

  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  
  UInt i_start = (UInt)(start/dx[0]);
  UInt i_end = (UInt)(end/dx[0]);

  UInt nb_dx = (UInt)((i_end-i_start)/(2*number));  
  UInt new_number = std::ceil((i_end-i_start)/(2*nb_dx));

  i_end = i_start+2*new_number*nb_dx;

  if(nb_dx < 2) {
    nb_dx = 1;
    std::cout << "! Asperity width was to small and is currently set to grid size !" << std::endl;
  }

  Real tough_ratio = (1.+delta_max_shr_strength)/(1.-delta_max_shr_strength);
  Real asperity_size = 2*nb_dx*dx[0];

  std::cout << "Asperity width (strong+weak) = " << asperity_size << std::endl;
  std::cout << "Asperity toughness ratio = " << tough_ratio << std::endl;
  std::cout << "Heterogeneous area extended up to " << i_end*dx[0] << std::endl;

  for (UInt ix = i_start; ix < i_end; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt no = ix-i_start;
      UInt asp_bool = ((no-no%nb_dx)%(2*nb_dx))/nb_dx;
      UInt index = ix+iz*n_ele[0];
      
      if ((asp_bool==1)) {

	(*nor_strength)[index] *= (1+delta_max_nor_strength);      
	(*shr_strength)[index] *= (1+delta_max_shr_strength);
	(*crit_n_open)[index] *= (1+delta_crit_nor_opening);
	(*crit_s_open)[index] *= (1+delta_crit_shr_opening);
	(*ind_crack)[index] = 5;
      }
      else {
	(*nor_strength)[index] *= (1-delta_max_nor_strength);      
	(*shr_strength)[index] *= (1-delta_max_shr_strength);
	(*crit_n_open)[index] *= (1-delta_crit_nor_opening);
	(*crit_s_open)[index] *= (1-delta_crit_shr_opening);
	(*ind_crack)[index] = 4;
      }
    }
  }

  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " HETEROGENEOUS INTERFACE " << std::endl
	      << "* Delta critical normal opening: " << delta_crit_nor_opening << std::endl
	      << "* Delta critical shear opening: " << delta_crit_shr_opening << std::endl
	      << "* Delta maximal normal strength: " << delta_max_nor_strength << std::endl
	      << "* Delta maximal shear strength: " << delta_max_shr_strength << std::endl
	      << std::endl;	

  out_parameters << "toughness_ratio " << tough_ratio << std::endl
		 << "asperity_paired_size " << asperity_size << std::endl
		 << "delta_delta_c_nor " << delta_crit_nor_opening << std::endl
		 << "delta_delta_c_shr " << delta_crit_shr_opening << std::endl
		 << "delta_tau_max_nor " << delta_max_nor_strength << std::endl
		 << "delta_tau_max_shr " << delta_max_shr_strength << std::endl;
  return i_end;
}

/* -------------------------------------------------------------------------- */
template<>
inline UInt Interfacer<_linear_coupled_cohesive>::createOneXStripe(Real start_x, Real end_x,
								   Real start_z, Real end_z,
								   Real delta_max_nor_strength, 
								   Real delta_max_shr_strength,
								   Real delta_crit_nor_opening, 
								   Real delta_crit_shr_opening) {

  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  
  UInt ix_start = (UInt)(start_x/dx[0]);
  UInt ix_end = (UInt)(end_x/dx[0]);
  
  UInt iz_start = (UInt)(start_z/dx[1]);
  UInt iz_end = (UInt)(end_z/dx[1]);
    
  if((iz_end-iz_start) < 1) {
    iz_end = iz_start + 1;
    std::cout << "! Asperity width in z was to small and is currently set to grid size !" << std::endl;
  }
  
  Real tough_ratio = (1.+abs(delta_max_shr_strength))/(1.-abs(delta_max_shr_strength));
  tough_ratio = tough_ratio * tough_ratio;
  Real asperity_size_z = (iz_end-iz_start)*dx[1];

  std::cout << "Asperity width = " << asperity_size_z << std::endl;
  std::cout << "This Asperity has a toughness ratio of = " << tough_ratio << std::endl;
  
  for (UInt ix = ix_start; ix < ix_end; ++ix) {
    for (UInt iz = iz_start; iz < iz_end; ++iz) {
      
      UInt index = ix+iz*n_ele[0];
      
      (*nor_strength)[index] *= (1+delta_max_nor_strength);      
      (*shr_strength)[index] *= (1+delta_max_shr_strength);
      (*crit_n_open)[index] *= (1+delta_crit_nor_opening);
      (*crit_s_open)[index] *= (1+delta_crit_shr_opening);
    }
  }
  
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " HETEROGENEOUS INTERFACE " << std::endl
	      << "* Delta critical normal opening: " << delta_crit_nor_opening << std::endl
	      << "* Delta critical shear opening: " << delta_crit_shr_opening << std::endl
	      << "* Delta maximal normal strength: " << delta_max_nor_strength << std::endl
	      << "* Delta maximal shear strength: " << delta_max_shr_strength << std::endl
	      << std::endl;
  
  //out_parameters << "toughness_ratio " << tough_ratio << std::endl
  out_parameters << "asperity_paired_size " << asperity_size_z << std::endl
		 << "delta_delta_c_nor " << delta_crit_nor_opening << std::endl
		 << "delta_delta_c_shr " << delta_crit_shr_opening << std::endl
		 << "delta_tau_max_nor " << delta_max_nor_strength << std::endl
		 << "delta_tau_max_shr " << delta_max_shr_strength << std::endl;
  return iz_end;
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughCenteredCrack(Real initial_crack_size,
									     Real crit_nor_opening, 
									     Real crit_shr_opening, 
									     Real max_nor_strength, 
									     Real max_shr_strength){
  
  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];
  
  std::fill(crit_n_open->begin(), crit_n_open->end(), crit_nor_opening);
  std::fill(crit_s_open->begin(), crit_s_open->end(), crit_shr_opening);
  std::fill(res_n_strength->begin(),  res_n_strength->end(), 0);
  std::fill(res_s_strength->begin(),  res_s_strength->end(), 0);
  
  Real pos = 0.5*dx[0];
   
  for (UInt ix = 0; ix < n_ele[0]/2; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt right_cmpnt = ix+n_ele[0]/2+iz*n_ele[0];
      UInt left_cmpnt = n_ele[0]/2-ix-1+iz*n_ele[0];

      if (pos <= initial_crack_size/2){

	(*nor_strength)[left_cmpnt]=0;
	(*shr_strength)[left_cmpnt]=0;
	(*ind_crack)[left_cmpnt] = 2;

	(*nor_strength)[right_cmpnt]=0;
	(*shr_strength)[right_cmpnt]=0;
	(*ind_crack)[right_cmpnt] = 2;
      }
      else {
	(*nor_strength)[left_cmpnt]=max_nor_strength;
	(*shr_strength)[left_cmpnt]=max_shr_strength;
	(*ind_crack)[left_cmpnt] = 0;

	(*nor_strength)[right_cmpnt]=max_nor_strength;
	(*shr_strength)[right_cmpnt]=max_shr_strength;
	(*ind_crack)[right_cmpnt] = 0;
      }
    }
    pos +=dx[0];
  }

  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " CENTERED CRACK " << std::endl
	      << "* Initial crack size: " << initial_crack_size << std::endl
	      << "* Critical normal opening: " << crit_nor_opening << std::endl
	      << "* Critical shear opening: " << crit_shr_opening << std::endl
	      << "* Maximal normal strength: " << max_nor_strength << std::endl
	      << "* Maximal shear strength: " << max_shr_strength << std::endl
	      << std::endl;	

  out_parameters << "a0 " << initial_crack_size << std::endl
		 << "delta_c_nor " << crit_nor_opening << std::endl
		 << "delta_c_shr " << crit_shr_opening << std::endl
		 << "tau_max_nor " << max_nor_strength << std::endl
		 << "tau_max_shr " << max_shr_strength << std::endl;  
}


/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughLeftSidedCrack(Real initial_crack_size,
									      Real crit_nor_opening, 
									      Real crit_shr_opening, 
									      Real max_nor_strength, 
									      Real max_shr_strength){
  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];

  std::fill(crit_n_open->begin(), crit_n_open->end(), crit_nor_opening);
  std::fill(crit_s_open->begin(), crit_s_open->end(), crit_shr_opening);
  std::fill(res_n_strength->begin(),  res_n_strength->end(), 0);
  std::fill(res_s_strength->begin(),  res_s_strength->end(), 0);
  
  Real pos = 0.5*dx[0];

  for (UInt ix = 0; ix < n_ele[0]; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt cmpnt = ix+iz*n_ele[0];
    
      if (pos <= initial_crack_size){
	(*nor_strength)[cmpnt]=0;
	(*shr_strength)[cmpnt]=0;
	(*ind_crack)[cmpnt] = 2;
      }
      else {
	(*nor_strength)[cmpnt]=max_nor_strength;
	(*shr_strength)[cmpnt]=max_shr_strength;
	(*ind_crack)[cmpnt] = 0;
      }
    }
    pos +=dx[0];
  }
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " LEFT-SIDED CRACK " << std::endl
	      << "* Initial crack size: " << initial_crack_size << std::endl
	      << "* Critical normal opening: " << crit_nor_opening << std::endl
	      << "* Critical shear opening: " << crit_shr_opening << std::endl
	      << "* Maximal normal strength: " << max_nor_strength << std::endl
	      << "* Maximal shear strength: " << max_shr_strength << std::endl
	      << std::endl;	
  
  out_parameters << "a0 " << initial_crack_size << std::endl
		 << "delta_c_nor " << crit_nor_opening << std::endl
		 << "delta_c_shr " << crit_shr_opening << std::endl
		 << "tau_max_nor " << max_nor_strength << std::endl
		 << "tau_max_shr " << max_shr_strength << std::endl;  
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createRightPropagatingCrackRoundAsp(Real initial_crack_size,
										      Real crit_nor_opening, 
										      Real crit_shr_opening,
										      Real max_nor_strength, 
										      Real max_shr_strength, 
										      Real radius, 
										      std::vector<Real> asp_ctr,
										      Real ratio_strength,
										      Real ratio_crit_open) {

  std::vector<Real> pos = {0.5*dx[0], 0.5*dx[1]}; 

  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];

  std::fill(crit_n_open->begin(), crit_n_open->end(), crit_nor_opening);
  std::fill(crit_s_open->begin(), crit_s_open->end(), crit_shr_opening);
  std::fill(res_n_strength->begin(),  res_n_strength->end(), 0);
  std::fill(res_s_strength->begin(),  res_s_strength->end(), 0);

  for (UInt i = 0; i < n_ele[0]; ++i) {
    pos[1]=0.5*dx[1];
    for (UInt j = 0; j < n_ele[1]; ++j) {
     
      UInt cmpnt = i+n_ele[0]*j;

      if (pos[0] <= initial_crack_size){
	(*nor_strength)[cmpnt]=0;
	(*shr_strength)[cmpnt]=0;
	(*ind_crack)[cmpnt] = 2;
      }
           
      //round Asperity
      else if ((pos[0]-asp_ctr[0])*(pos[0]-asp_ctr[0])+(pos[1]-asp_ctr[1])*(pos[1]-asp_ctr[1])<radius*radius){

	(*nor_strength)[cmpnt]=ratio_strength*max_nor_strength;
	(*shr_strength)[cmpnt]=ratio_strength*max_shr_strength;
	
	(*crit_n_open)[cmpnt] = ratio_crit_open*crit_nor_opening;
	(*crit_s_open)[cmpnt] = ratio_crit_open*crit_shr_opening;

	(*ind_crack)[cmpnt]=5;
	}
      
      else{
	
	(*nor_strength)[cmpnt]=max_nor_strength;
	(*shr_strength)[cmpnt]=max_shr_strength;
	(*ind_crack)[cmpnt] = 0;
      }
   
      pos[1] += dx[1];
    }
    pos[0] += dx[0];
  }

  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " CIRCULAR ASPERITY " << std::endl
	      << "* Initial crack size: " << initial_crack_size << std::endl
	      << "* Critical normal opening: " << crit_nor_opening << std::endl
	      << "* Critical shear opening: " << crit_shr_opening << std::endl
	      << "* Maximal normal strength: " << max_nor_strength << std::endl
	      << "* Maximal shear strength: " << max_shr_strength << std::endl
	      << "* Asperity center (x,z): (" << asp_ctr[0] << "," << asp_ctr[1] << ")" << std::endl
	      << "* Asperity radius: " << radius << std::endl
	      << "* Strength ratio: " << ratio_strength << std::endl
    	      << "* Critical opening ratio: " << ratio_crit_open << std::endl
	      << std::endl;	

  out_parameters << "a0 " << initial_crack_size << std::endl
		 << "delta_c_nor " << crit_nor_opening << std::endl
		 << "delta_c_shr " << crit_shr_opening << std::endl
		 << "tau_max_nor " << max_nor_strength << std::endl
		 << "tau_max_shr " << max_shr_strength << std::endl
    		 << "ratio_str " << ratio_strength << std::endl
    		 << "ratio_delta_c " << ratio_crit_open << std::endl;
}

/* -------------------------------------------------------------------------- */
template<>
void Interfacer<_linear_coupled_cohesive>::createIncohIntfc() {

  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * res_n_strength = datas[_residual_normal_strength];
  std::vector<Real> * res_s_strength = datas[_residual_shear_strength];

  std::fill(nor_strength->begin(),  nor_strength->end(), 0);
  std::fill(shr_strength->begin(),  shr_strength->end(), 0);
  std::fill(ind_crack->begin(),  ind_crack->end(), 2);
  std::fill(res_n_strength->begin(),  res_n_strength->end(), 0);
  std::fill(res_s_strength->begin(),  res_s_strength->end(), 0);
  
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " INTERFACE WITHOUT INITIAL COHESION" << std::endl;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* ------------------TEMPLATES FOR VISCOELASTIC INTERFACES       ------------ */
/* ------------------SHOULD BE REPLACED BY A GENERAL TEMPLATE WITH SWITCH---- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_viscoelastic_coupled_cohesive>::createUniformInterface() {

  Real crit_nor_opening = DataRegister::getParameter<Real>("critical_normal_opening");
  Real crit_shr_opening = DataRegister::getParameter<Real>("critical_shear_opening");
  Real max_nor_strength = DataRegister::getParameter<Real>("max_normal_strength");
  Real max_shr_strength = DataRegister::getParameter<Real>("max_shear_strength");
  Real lim_velocity = DataRegister::getParameter<Real>("lim_velocity");

  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * lim_vel = datas[_lim_velocity];
  
  std::fill(nor_strength->begin(), nor_strength->end(), max_nor_strength);
  std::fill(shr_strength->begin(), shr_strength->end(), max_shr_strength);
  std::fill(crit_n_open->begin(), crit_n_open->end(), crit_nor_opening);
  std::fill(crit_s_open->begin(), crit_s_open->end(), crit_shr_opening);
  std::fill(lim_vel->begin(), lim_vel->end(), lim_velocity);
  
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " VISCOELASTIC COHESIVE LAW " << std::endl
	      << " UNIFORM INTERFACE " << std::endl
	      << "* Critical normal opening: " << crit_nor_opening << std::endl
	      << "* Critical shear opening: " << crit_shr_opening << std::endl
	      << "* Maximal normal strength: " << max_nor_strength << std::endl
	      << "* Maximal shear strength: " << max_shr_strength << std::endl
	      << "* Limit velocity: " << lim_velocity << std::endl
	      << std::endl;	

  out_parameters << "delta_c_nor " << crit_nor_opening << std::endl
		 << "delta_c_shr " << crit_shr_opening << std::endl
		 << "tau_max_nor " << max_nor_strength << std::endl
		 << "tau_max_shr " << max_shr_strength << std::endl
		 << "lim_vel " << lim_velocity << std::endl;  
}
/* -------------------------------------------------------------------------- */

template<>
inline void Interfacer<_viscoelastic_coupled_cohesive>::createThroughArea(Real area_start, Real area_end,
									  UInt cracking_index,
									  Real ratio_max_nor_strength, 
									  Real ratio_max_shr_strength,
									  Real ratio_crit_nor_opening, 
									  Real ratio_crit_shr_opening,
									  bool vrtr) {
  
  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];

  UInt i_start = (UInt)(area_start/dx[0]);
  UInt i_end = (UInt)(area_end/dx[0]);

  for (UInt ix = i_start; ix < i_end; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt index = ix+iz*n_ele[0];

      (*nor_strength)[index] *= (vrtr+ratio_max_nor_strength);      
      (*shr_strength)[index] *= (vrtr+ratio_max_shr_strength);
      (*crit_n_open)[index] *= (vrtr+ratio_crit_nor_opening);
      (*crit_s_open)[index] *= (vrtr+ratio_crit_shr_opening);
      (*ind_crack)[index] = cracking_index;
    }
  }
}

/* -------------------------------------------------------------------------- */

template<>
inline void Interfacer<_viscoelastic_coupled_cohesive>::createThroughCrack(Real crack_start, Real crack_end) {
  
  createThroughArea(crack_start, crack_end, 2, 0., 0.);
  Real initial_crack_size = crack_end-crack_start;
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " THROUGH CRACK " << std::endl
	      << "* Initial crack size: " << initial_crack_size << std::endl
	      << "* Crack starts at: " << crack_start << std::endl
	      << "* Crack ends at: " << crack_end << std::endl
	      << std::endl;	

  out_parameters << "a0 " << initial_crack_size << std::endl
		 << "a_l " << crack_start << std::endl
		 << "a_r " << crack_end << std::endl;
}


/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_viscoelastic_coupled_cohesive>::createThroughWall(Real wall_start, Real wall_end) {
  
  createThroughArea(wall_start, wall_end, 6, 1000., 1000.);

  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " THROUGH WALL " << std::endl
	      << "* Wall starts at: " << wall_start << std::endl
	      << "* Wall ends at: " << wall_end << std::endl
	      << std::endl;	
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_viscoelastic_coupled_cohesive>::createHeterogeneousInterface(std::vector<Real> crit_nor_opening, std::vector<Real> max_nor_strength,std::vector<Real> res_nor_strength, std::vector<Real> crit_shr_opening, std::vector<Real> max_shr_strength, std::vector<Real> res_shr_strength) {

  Real lim_velocity = DataRegister::getParameter<Real>("lim_velocity");
  
  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  std::vector<Real> * crit_n_open = datas[_critical_normal_opening];
  std::vector<Real> * crit_s_open = datas[_critical_shear_opening];
  std::vector<Real> * lim_vel = datas[_lim_velocity];

  // Fill the array containing the limiting velocity
  std::fill(lim_vel->begin(), lim_vel->end(), lim_velocity);
  
  // Check the length of the shear arguments...
  if ( crit_shr_opening.size() == 0 ){
    crit_shr_opening = crit_nor_opening;
    max_shr_strength = max_nor_strength;
  }
  
  for (UInt x = 0; x < n_ele[0]; ++x) {
    for (UInt z = 0; z < n_ele[1]; ++z ) {
      
      (*shr_strength)[x+z*n_ele[0]] = max_shr_strength[x+z*n_ele[0]];
      (*nor_strength)[x+z*n_ele[0]] = max_nor_strength[x+z*n_ele[0]];
      (*crit_s_open)[x+z*n_ele[0]] = crit_shr_opening[x+z*n_ele[0]];
      (*crit_n_open)[x+z*n_ele[0]] = crit_nor_opening[x+z*n_ele[0]];
      
      if(max_shr_strength[x+z*n_ele[0]]==0 || max_nor_strength[x+z*n_ele[0]]) {
	(*ind_crack)[x+z*n_ele[0]] = 0;
      }
    }
  }
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << "Maximum strength and critical opening patterns inserted from vectors."
	      << std::endl
	      << std::endl;
}
