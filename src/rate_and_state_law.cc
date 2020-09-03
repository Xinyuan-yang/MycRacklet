/* -------------------------------------------------------------------------- */
#include "rate_and_state_law.hh"
#include <fstream>
#include <algorithm>
#include <random>
#include <iomanip>
/* -------------------------------------------------------------------------- */
void RateAndStateLaw::initStateEvolution() {
  state_evol = std::make_shared<StateEvolution>();
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::initRegularizedStateEvolution(Real v0) {
  state_evol = std::make_shared<RegularizedStateEvolution>(v0);
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
void RateAndStateLaw::initRegularizedFormulation(Real v0, Real theta, Real xi) {
  formulation = std::make_shared<RegularizedRandSFormulation>(v0,theta,xi);
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::initRegularizedWeakeningFormulation(Real v0, Real theta, Real xi) {
  formulation = std::make_shared<RegularizedWeakeningRandSFormulation>(v0,theta,xi);
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::setVelocityPredictor(std::vector<Real> v_0_pred) {

  for (UInt h = 0; h < n_ele[0]; ++h) {
    for (UInt j = 0; j < n_ele[1]; ++j) {
      for (UInt dir=0; dir < dim; ++dir) {
	UInt i=h+j*n_ele[0];
	(*dot_u_top)[i*dim+dir] = 0.5*v_0_pred[dir]/c_s;
	(*dot_u_bot)[i*dim+dir] = -0.5*v_0_pred[dir]/c_s;
      }
    }
  }
}
/* -------------------------------------------------------------------------- */
void RateAndStateLaw::setV0(std::vector<Real> v_0) {

  std::cout << "SETTING V0" << std::endl;

  for (UInt dir=0; dir < dim; ++dir) {
    V_0[dir] = v_0[dir];
  }

  Real V_norm = sqrt(V_0[0]*V_0[0]+V_0[1]*V_0[1]);
  std::cout << "Steady-state initial velocity V="<< V_norm<< std::endl;
  DataRegister::out_parameters << std::setprecision(15) << "V0 " << V_norm << std::endl;

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
    Real rate = (*shear_velo_jump)[i]*c_s;
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
      V_0[s] = rate*shear_ratio[s];
      (*intfc_trac)[3*i+2*s] = stress*shear_ratio[s];
    }
    
    fric_strength[i] = sigma_0*cf[i];
  }
  Real V_norm = sqrt(V_0[0]*V_0[0]+V_0[1]*V_0[1]);
  std::cout << "Steady-state initial velocity V="<<V_norm<<" [m/s] with friction coefficient cf="
	    << cf[0] << std::endl;
  DataRegister::out_parameters << "V0 " << V_norm << std::endl;
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

  Real V0 = sqrt(V_0[0]*V_0[0]+V_0[1]*V_0[1]);
  
  for (UInt i = 0; i < n_ele; ++i) {     

    Real delta_cf,K;
    Real gamma = 0.1;
    UInt nb_t_stp = 1;
    bool is_good = false;
    Real stress, rate;
    std::vector<Real> shear_ratio = {tau_x[i]/shear_stress[i],tau_z[i]/shear_stress[i]};
    
    while(!is_good) {
    
      stress = shear_stress[i];    
      rate = (*shear_velo_jump)[i]*c_s + V0;
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
	//phi[i]=new_phi;
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
      av_vel[s] += 0.5*(rate*shear_ratio[s]-V_0[s]);
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

  Real V0 = sqrt(V_0[0]*V_0[0]+V_0[1]*V_0[1]);
  
  for (UInt i = 0; i < n_ele; ++i) {     

    Real delta_cf,K;
    Real gamma = 0.1;
    UInt nb_t_stp = 1;
    bool is_good = false;
    Real stress, rate;
    std::vector<Real> shear_ratio = {tau_x[i]/shear_stress[i],tau_z[i]/shear_stress[i]};
    
    while(!is_good) {
    
      stress = shear_stress[i];    
      rate = (*shear_velo_jump)[i]*c_s + V0;
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
      (*dot_u_top)[i*3+2*s] = 0.5*(rate*shear_ratio[s]-V_0[s])/c_s;
      (*dot_u_bot)[i*3+2*s] = -0.5*(rate*shear_ratio[s]-V_0[s])/c_s;
      (*intfc_trac)[3*i+2*s] = stress*shear_ratio[s] - accoust*(rate*shear_ratio[s]-V_0[s])/2;
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
void RateAndStateLaw::insertPerturbationPatch(std::vector<UInt> patch_limits, Real new_rate) {
  
  UInt n_ele = D.size();
  CrackProfile * stresses = datas[_top_dynamic_stress];
  
  if(patch_limits[1]>(n_ele-1))
    cRacklet::error("Perturbation patch exceeds the size of the interface");
  
  for (UInt i = patch_limits[0]; i < patch_limits[1]; ++i) {     

    (*dot_u_top)[i*3+2] = 0.5*(new_rate-V_0[1])/c_s;
    (*dot_u_bot)[i*3+2] = -0.5*(new_rate-V_0[1])/c_s;

    Real tractions = (*stresses)[i*3+2] - accoust*(new_rate-V_0[1])/2;
    
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
    Real new_rate = amplitude*gauss+V_0[1];

    (*dot_u_top)[i*3+2] = 0.5*(new_rate-V_0[1])/c_s;
    (*dot_u_bot)[i*3+2] = -0.5*(new_rate-V_0[1])/c_s;

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
  
}

/* -------------------------------------------------------------------------- */
