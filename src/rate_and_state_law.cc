/**
 * @file   rate_and_state_law.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Fri Feb 24 16:41:05 2017
 *
 * @brief  Implementation of the RateAndStateLaw class
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
#include "rate_and_state_law.hh"
#include <fstream>
#include <algorithm>
#include <random>
#include <iomanip>
#include <math.h>
/* -------------------------------------------------------------------------- */
void RateAndStateLaw::initStateEvolution() {
  state_evol = std::make_shared<StateEvolution>();
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::initRegularizedStateEvolution(Real v0) {
  state_evol = std::make_shared<RegularizedStateEvolution>(v0);
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::initSlipStateEvolution() {
  state_evol = std::make_shared<SlipStateEvolution>();
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::initStandardFormulation() {
  formulation = std::make_shared<RandSFormulation>();
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::initVelocityWeakeningFormulation() {
  formulation = std::make_shared<WeakeningRandSFormulation>();
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::initRegularizedFormulation(Real v0) {
  formulation = std::make_shared<RegularizedRandSFormulation>(v0);
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::initRegularizedWeakeningFormulation(Real v0) {
  formulation = std::make_shared<RegularizedWeakeningRandSFormulation>(v0);
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::setVelocityPredictor(std::vector<Real> v_0_pred) {

  for (UInt h = 0; h < n_ele[0]; ++h) {
    for (UInt j = 0; j < n_ele[1]; ++j) {
      for (UInt dir=0; dir < dim; ++dir) {
	UInt i=h+j*n_ele[0];
	(*dot_u_top)[i*dim+dir] = 0.5*v_0_pred[dir];
	(*dot_u_bot)[i*dim+dir] = -0.5*v_0_pred[dir];
      }
    }
  }
}
/* -------------------------------------------------------------------------- */
void RateAndStateLaw::setV0(std::vector<Real> v_0) {

  std::cout << "SETTING V0" << std::endl;

  for (UInt h = 0; h < n_ele[0]; ++h) {
    for (UInt j = 0; j < n_ele[1]; ++j) {
      for (UInt dir=0; dir < dim; ++dir) {
      	UInt i=h+j*n_ele[0];
	V_0[i*2+dir] = v_0[dir];
      }
    }
  }

  Real V_norm = sqrt(V_0[0]*V_0[0]+V_0[1]*V_0[1]);
  std::cout << "Steady-state initial velocity of the first element V="<< V_norm<< std::endl;
  DataRegister::out_parameters << std::setprecision(15) << "V0 " << V_norm << std::endl;

  saveV0();
  
}
/* -------------------------------------------------------------------------- */
void RateAndStateLaw::saveV0() {

  std::cout << "Saving V0 to a file"<<std::endl;
  
  std::string output_v0 = output_dir+"v0_profile.cra";
  std::ofstream outFile(output_v0);
  for (std::vector<Real>::iterator it = V_0.begin(); it != V_0.end(); ++it){
    outFile << *it << "\n";
  }
}
/* -------------------------------------------------------------------------- */
void RateAndStateLaw::initInterfaceConditions() {

  computeSteadyStateSliding();
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::computeSteadyStateSliding() {
  
  UInt n_ele = D.size();
  
  CrackProfile tau_x = stresses->getStridedPart(3,0);
  CrackProfile tau_z = stresses->getStridedPart(3,2);

  CrackProfile shear_stress = tau_x*tau_x + tau_z*tau_z;
  shear_stress.squareRoot();

  for (UInt i = 0; i < n_ele; ++i) {     

    Real stress = shear_stress[i];    
    Real rate = (*shear_velo_jump)[i];
    Real converg=1.0;
    Real delta_v = 0.0;
    Real delta_cf,K;
    Real n=0;
    std::vector<Real> shear_ratio = {tau_x[i]/stress,tau_z[i]/stress};
    
    while ((!cRacklet::has_converged(converg))) {

      rate += delta_v;
      
      if(rate == 0.)
	phi[i] = 0.;
      else
	phi[i] = state_evol->getSteadyState(rate,D[i]);
      cf[i] = (*formulation)(rate,phi[i],a[i],b[i],D[i],f_0[i],v_star[i],phi_star[i]);
      delta_cf = stress/sigma_0 - cf[i];
      K = formulation->getSteadyTangent(rate,phi[i],a[i],b[i],D[i],f_0[i],v_star[i],phi_star[i]);
      
      delta_v = delta_cf/K;
      if ((delta_v > 0))
	delta_v = std::min(std::abs(delta_v),0.1*std::abs(rate));
      else
	delta_v = -std::min(std::abs(delta_v),0.1*std::abs(rate));
      converg = delta_v/rate;
      n += 1;
    }
    
    for (UInt s = 0; s < 2; ++s) {
      (*dot_u_top)[i*3+2*s] = 0.;
      (*dot_u_bot)[i*3+2*s] = 0.;
      (*intfc_trac)[3*i+2*s] = stress*shear_ratio[s];

      V_0[i*2+s] = rate*shear_ratio[s];
    }
    
    fric_strength[i] = sigma_0*cf[i];
  }
  Real V_norm = sqrt(V_0[0]*V_0[0]+V_0[1]*V_0[1]);
  std::cout << "Steady-state initial velocity of the first element V="<< V_norm<<" [m/s] with friction coefficient cf="
	    << cf[0] << std::endl;
  DataRegister::out_parameters << "V0 " << V_norm << std::endl;

  saveV0();
  
}

/* -------------------------------------------------------------------------- */
std::vector<Real> RateAndStateLaw::computeNextAverageVelocity() {
  
  std::vector<Real> av_vel;
  av_vel.resize(2);
  
  UInt n_ele = D.size();
  
  CrackProfile tau_x = stresses->getStridedPart(3,0);
  CrackProfile tau_z = stresses->getStridedPart(3,2);
  
  CrackProfile shear_stress = tau_x*tau_x + tau_z*tau_z;
  shear_stress.squareRoot();

  
  for (UInt i = 0; i < n_ele; ++i) {     

    Real V0 = sqrt(V_0[i*2+0]*V_0[i*2+0]+V_0[i*2+1]*V_0[i*2+1]);
    
    Real delta_cf,K;
    Real gamma = 0.1;
    UInt nb_t_stp = 1;
    bool is_good = false;
    Real stress, rate;
    std::vector<Real> shear_ratio = {tau_x[i]/shear_stress[i],tau_z[i]/shear_stress[i]};
    
    while(!is_good) {
      
      stress = shear_stress[i];    
      rate = (*shear_velo_jump)[i] + V0;
      Real new_phi = phi[i];
      
      for (UInt t = 0; t < nb_t_stp; ++t) {

	Real converg_v = 1.0;
	Real converg_cf = 1.0;
	Real delta_v = 0.0;
	UInt n=0;

	while ((!cRacklet::has_converged(converg_v))&&(!cRacklet::has_converged(converg_cf))) { //&&(is_good)) {
	
	  rate += delta_v;
     
	  cf[i] = (*formulation)(rate,new_phi,a[i],b[i],D[i],f_0[i],v_star[i],phi_star[i]); 
	  delta_cf = stress/sigma_0 - cf[i] - accoust*(rate-V0)/(2*sigma_0);
	  K = accoust/(2*sigma_0) + formulation->getTangent(rate,new_phi,a[i],b[i],D[i],f_0[i],v_star[i],phi_star[i]);
	  delta_v = delta_cf/K;
	  if ((delta_v > 0))
	    delta_v = std::min(std::abs(delta_v),gamma*std::abs(rate));
	  else
	    delta_v = -std::min(std::abs(delta_v),gamma*std::abs(rate));
	  converg_v = delta_v/rate;
	  converg_cf = delta_cf;
	  n += 1;
	}
	new_phi += (*state_evol)(rate, new_phi, D[i], delta_t/(Real)(nb_t_stp));
      }
      
      if(new_phi<0) {
	nb_t_stp *= 2;
	is_good = false;
      }
      else {
	is_good = true;
      }
      
      if(nb_t_stp>1) {
	std::cerr << "WARNING ! Negative Phi values have been observed ! " 
		  << "Internal time step is locally refined, leading to potential instabilities. " 
		  << std::endl;
      }
    }
    
    if(phi[i]<0) {
      std::cout << "Phi: " << phi[i] << " V: " << rate
		<< std::endl;
      abort();
    } 
    
    for (UInt s = 0; s < 2; ++s) {
      av_vel[s] += 0.5*(rate*shear_ratio[s]-V_0[i*2+s]);
    }
  }      
  
  for (UInt s=0; s < 2; ++s){
    av_vel[s] /= n_ele;
  }
  
  return av_vel;
  
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::updateInterfaceConditions() {
  
  UInt n_ele = D.size();
  
  CrackProfile tau_x = stresses->getStridedPart(3,0);
  CrackProfile tau_z = stresses->getStridedPart(3,2);
  
  CrackProfile shear_stress = tau_x*tau_x + tau_z*tau_z;
  shear_stress.squareRoot();

  
  for (UInt i = 0; i < n_ele; ++i) {     

    Real V0 = sqrt(V_0[i*2+0]*V_0[i*2+0]+V_0[i*2+1]*V_0[i*2+1]);
    
    Real delta_cf,K;
    Real gamma = 0.1;
    UInt nb_t_stp = 1;
    bool is_good = false;
    Real stress, rate;
    std::vector<Real> shear_ratio = {tau_x[i]/shear_stress[i],tau_z[i]/shear_stress[i]};
    
    while(!is_good) {
    
      stress = shear_stress[i];    
      rate = (*shear_velo_jump)[i] + V0;
      Real new_phi = phi[i];
      
      for (UInt t = 0; t < nb_t_stp; ++t) {

	Real converg_v = 1.0;
	Real converg_cf = 1.0;
	Real delta_v = 0.0;
	UInt n=0;

	while ((!cRacklet::has_converged(converg_v))&&(!cRacklet::has_converged(converg_cf))) { //&&(is_good)) {
	
	  rate += delta_v;
     
	  cf[i] = (*formulation)(rate,new_phi,a[i],b[i],D[i],f_0[i],v_star[i],phi_star[i]); 
	  delta_cf = stress/sigma_0 - cf[i] - accoust*(rate-V0)/(2*sigma_0);
	  K = accoust/(2*sigma_0) + formulation->getTangent(rate,new_phi,a[i],b[i],D[i],f_0[i],v_star[i],phi_star[i]);
	  delta_v = delta_cf/K;
	  if ((delta_v > 0))
	    delta_v = std::min(std::abs(delta_v),gamma*std::abs(rate));
	  else
	    delta_v = -std::min(std::abs(delta_v),gamma*std::abs(rate));
	  converg_v = delta_v/rate;
	  converg_cf = delta_cf;
	  n += 1;
	}
	new_phi += (*state_evol)(rate, new_phi, D[i], delta_t/(Real)(nb_t_stp));
      }
           
      if(new_phi<0) {
	nb_t_stp *= 2;
	is_good = false;
      }
      else {
	is_good = true;
	phi[i]=new_phi;
      }

      if(nb_t_stp>1) {
	std::cerr << "WARNING ! Negative Phi values have been observed ! " 
		  << "Internal time step is locally refined, leading to potential instabilities. " 
		  << std::endl;
      }
    }
    
    if(phi[i]<0) {
      std::cout << "Phi: " << phi[i] << " V: " << rate
		<< std::endl;
      abort();
    } 

    for (UInt s = 0; s < 2; ++s) {
      (*dot_u_top)[i*3+2*s] = 0.5*(rate*shear_ratio[s]-V_0[i*2+s]);
      (*dot_u_bot)[i*3+2*s] = -0.5*(rate*shear_ratio[s]-V_0[i*2+s]);
      (*intfc_trac)[3*i+2*s] = stress*shear_ratio[s] - accoust*(rate*shear_ratio[s]-V_0[i*2+s])/2;
    }

    fric_strength[i] = sigma_0*cf[i];
  }      
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::perturbState(Real epsilon, Real k) {
  
  UInt n_ele = D.size();

  for (UInt i = 0; i < n_ele; ++i) {     
    phi[i] += epsilon*sin(k*i*2*M_PI/(Real)(n_ele));
  }
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::perturbState(std::vector<Real> perturbation) {
  
  UInt n_ele = D.size();

  for (UInt i = 0; i < n_ele; ++i) {     
    phi[i] += perturbation[i];
  }
}
/* -------------------------------------------------------------------------- */
void RateAndStateLaw::setStableState() {

  UInt n_ele = D.size();
  CrackProfile * stresses = datas[_top_dynamic_stress];

  for (UInt i = 0; i < n_ele; ++i) {

    Real slip_rate = 2*(*dot_u_top)[i*3+2]+V_0[i*2+1];
    Real tractions = (*stresses)[i*3+2] - accoust*(*dot_u_top)[i*3+2];

    phi[i] = formulation->getStableState(tractions,sigma_0,slip_rate,a[i],b[i],D[i],f_0[i],v_star[i],phi_star[i]);
    if (phi[0]<0)
      cRacklet::error("The current amplitude of the sliding velocity implies negative value of phi !");
  }
}
/* -------------------------------------------------------------------------- */
void RateAndStateLaw::insertPerturbationPatch(std::vector<UInt> patch_limits, Real new_rate) {
  
  UInt n_ele = D.size();
  CrackProfile * stresses = datas[_top_dynamic_stress];
  
  if(patch_limits[1]>(n_ele-1))
    cRacklet::error("Perturbation patch exceeds the size of the interface");
  
  for (UInt i = patch_limits[0]; i < patch_limits[1]; ++i) {     

    (*dot_u_top)[i*3+2] = 0.5*(new_rate-V_0[i*2+1]);
    (*dot_u_bot)[i*3+2] = -0.5*(new_rate-V_0[i*2+1]);

    Real tractions = (*stresses)[i*3+2] - accoust*(new_rate-V_0[i*2+1])/2;
    
    phi[i] = formulation->getStableState(tractions,sigma_0,new_rate,a[i],b[i],D[i],f_0[i],v_star[i],phi_star[i]);
    if (phi[0]<0)
      cRacklet::error("The amplitude of the perturbation patch implies negative value of phi !");
  }
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::insertGaussianPerturbation(Real std_dev, Real amplitude) {
  
  UInt n_ele = D.size();
  Real mu = n_ele/2.;
  CrackProfile * stresses = datas[_top_dynamic_stress];
  
  for (UInt i = 0; i < n_ele; ++i) {     

    Real gauss = exp(-0.5*((i-mu)/std_dev)*((i-mu)/std_dev));
    Real new_rate = amplitude*gauss+V_0[i*2+1];

    (*dot_u_top)[i*3+2] = 0.5*(new_rate-V_0[i*2+1]);
    (*dot_u_bot)[i*3+2] = -0.5*(new_rate-V_0[i*2+1]);

    Real tractions = (*stresses)[i*3+2] - accoust*(new_rate-V_0[i*2+1])/2;
    
    phi[i] = formulation->getStableState(tractions,sigma_0,new_rate,a[i],b[i],D[i],f_0[i],v_star[i],phi_star[i]);
    if (phi[0]<0)
      cRacklet::error("The amplitude of the gaussian perturbation implies negative value of phi !");
  }
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::insertSkewedPerturbation(Real std_dev, Real amplitude, Real alpha, Real rel_loc) {
  
  UInt n_ele = D.size();
  Real mu = n_ele*rel_loc;
  CrackProfile * stresses = datas[_top_dynamic_stress];

  std::vector<Real> skew_perturbation(n_ele); 
  
  for (UInt i = 0; i < n_ele; ++i) {     
    Real normal_pdf = exp(-0.5*((i-mu)/std_dev)*((i-mu)/std_dev));
    Real skew_cdf = 0.5*(1 + erf(alpha*(i-mu)/std_dev/sqrt(2)));
    skew_perturbation[i] = normal_pdf*skew_cdf;
  }

  // Compute max of the distribution

  Real max_skew = *std::max_element(skew_perturbation.begin(),skew_perturbation.end());
  
  // Second loop, add the perturbation
  for (UInt i = 0; i < n_ele; ++i) {     
    
    Real new_rate = amplitude*skew_perturbation[i]/max_skew+V_0[1];
    
    (*dot_u_top)[i*3+2] = 0.5*(new_rate-V_0[1]);
    (*dot_u_bot)[i*3+2] = -0.5*(new_rate-V_0[1]);
    
    Real tractions = (*stresses)[i*3+2] - accoust*(new_rate-V_0[1])/2;
    
    phi[i] = formulation->getStableState(tractions,sigma_0,new_rate,a[i],b[i],D[i],f_0[i],v_star[i],phi_star[i]);
    if (phi[0]<0)
      cRacklet::error("The amplitude of the gaussian perturbation implies negative value of phi !");
  }
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::insertSkewedPerturbation(std::vector<Real> std_dev, std::vector<Real> amplitude, std::vector<Real> alpha, std::vector<Real> rel_loc) {

  // This function adds several perturbations on top of each other
  
  UInt n_ele = D.size();
  std::vector<Real> mu = rel_loc;
  std::transform(mu.begin(), mu.end(), mu.begin(), [n_ele](Real &c){ return c*n_ele; });
  CrackProfile * stresses = datas[_top_dynamic_stress];

  std::vector<Real> total_skew_perturbation(n_ele); 
  
  for (UInt j = 0; j < std_dev.size(); ++j){
    std::vector<Real> skew_perturbation(n_ele); 
    for (UInt i = 0; i < n_ele; ++i) {     
      Real normal_pdf = exp(-0.5*((i-mu[j])/std_dev[j])*((i-mu[j])/std_dev[j]));
      Real skew_cdf = 0.5*(1 + erf(alpha[j]*(i-mu[j])/std_dev[j]/sqrt(2)));
      skew_perturbation[i] = normal_pdf*skew_cdf;
    }
    Real max_skew = *std::max_element(skew_perturbation.begin(),skew_perturbation.end());

    // Normalize the perturbation
    std::transform(skew_perturbation.begin(), skew_perturbation.end(), skew_perturbation.begin(), [j,amplitude,max_skew](Real &c){ return c*amplitude[j]/max_skew;});

    // Add to the total profile
    std::transform(total_skew_perturbation.begin(),total_skew_perturbation.end(),skew_perturbation.begin(),total_skew_perturbation.begin(),std::plus<Real>());
  }

  // Second loop, add the perturbation
  for (UInt i = 0; i < n_ele; ++i) {     
    
    Real new_rate = total_skew_perturbation[i]+V_0[1];
    
    (*dot_u_top)[i*3+2] = 0.5*(new_rate-V_0[1]);
    (*dot_u_bot)[i*3+2] = -0.5*(new_rate-V_0[1]);
    
    Real tractions = (*stresses)[i*3+2] - accoust*(new_rate-V_0[1])/2;
    
    phi[i] = formulation->getStableState(tractions,sigma_0,new_rate,a[i],b[i],D[i],f_0[i],v_star[i],phi_star[i]);
    if (phi[0]<0)
      cRacklet::error("The amplitude of the gaussian perturbation implies negative value of phi !");
  }
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::addGaussianNoiseToStateField(Real std_dev) {
  
  UInt n_ele = D.size();

  Real mean = 0.0;
  std::default_random_engine generator;
  std::normal_distribution<Real> dist(mean,std_dev);

  for (UInt i = 0; i < n_ele; ++i) {     
    phi[i] += dist(generator)*phi[i];
  }
  
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::insertPerturbationFromFile(std::string input_file) {
  
  UInt n_ele = D.size();

  std::ifstream fin(input_file); //opening the perturbation file

  Real pert;

  for (UInt i = 0; i < n_ele; ++i) {     
    fin >> pert;
    phi[i] += pert*phi[i];
  }
  
}

/* -------------------------------------------------------------------------- */

void RateAndStateLaw::restart(bool pausing, UInt nele_2d) {
  
  std::vector<Real> * state = datas[_state_variable];
  DataRegister::restartData(*state,"restart_state.cra",pausing, nele_2d);
  DataRegister::restartData(V_0,"restart_V0.cra",pausing, nele_2d);
}

/* -------------------------------------------------------------------------- */
