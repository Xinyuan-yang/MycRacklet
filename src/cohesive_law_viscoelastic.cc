/* -------------------------------------------------------------------------- */
#include "cohesive_law_viscoelastic.hh"
#include <fstream>
#include <algorithm>
#include <math.h> 
/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::initLinearFormulation() {
  formulation = std::make_shared<ViscoelasticFormulation>();
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::initPowerLawFormulation() {
  formulation = std::make_shared<ViscoelasticPowerLawFormulation>();
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::initInterfaceConditions() {
  
  memcpy(&this->max_nor_strength[0],&this->nor_strength[0],this->nor_strength.size()*sizeof(Real));
  memcpy(&this->max_shr_strength[0],&this->shr_strength[0],this->shr_strength.size()*sizeof(Real));
  memcpy(&this->lim_velocity[0],&this->lim_velocity[0],this->lim_velocity.size()*sizeof(Real));
  
  this->computeInitialVelocities();
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::restart(bool pausing, UInt nele_2d) {

  std::vector<Real> * nor_strength = datas[_normal_strength];
  DataRegister::restartData(*nor_strength,"restart_normal_strength.cra",pausing, nele_2d);
  std::vector<Real> * shr_strength = datas[_shear_strength];
  DataRegister::restartData(*shr_strength,"restart_shear_strength.cra",pausing, nele_2d);
  
  if(contact_law)
    contact_law->restart(pausing,nele_2d);
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::computeInitialVelocities() {

  Real strength;
  Real shr_trac;
  Real shr_velo;
  
  // Get the stresses

  std::vector<CrackProfile*> loads = {datas[_top_loading],datas[_bottom_loading]};
  
  CrackProfile load_x_top = loads[0]->getStridedPart(3,0);
  CrackProfile load_y_top = loads[0]->getStridedPart(3,1);
  CrackProfile load_z_top = loads[0]->getStridedPart(3,2);

  CrackProfile tau_x = load_x_top;
  CrackProfile tau_y = load_y_top;
  CrackProfile tau_z = load_z_top;

  
  CrackProfile v_x_top = velocities[0]->getStridedPart(3,0);
  CrackProfile v_y_top = velocities[0]->getStridedPart(3,1);
  CrackProfile v_z_top = velocities[0]->getStridedPart(3,2);

  CrackProfile v_x_bot = velocities[1]->getStridedPart(3,0);
  CrackProfile v_y_bot = velocities[1]->getStridedPart(3,1);
  CrackProfile v_z_bot = velocities[1]->getStridedPart(3,2);
    
  CrackProfile shear_stress = tau_x*tau_x + tau_z*tau_z;
  shear_stress.squareRoot();
  
  CrackProfile total_stress = tau_y * tau_y + shear_stress * shear_stress;;
  total_stress.squareRoot();
  
  UInt n=0;
  for (UInt x = 0; x < n_ele[0]; ++x) {
    for (UInt z = 0; z < n_ele[1]; ++z) {
      n=x+z*n_ele[0];

      bool broken_nor = false;  
      bool broken_shr = false;  
      
      // Check if the element is broken AND in compression.
      if((this->nor_strength[n]==0)&&(tau_y[n] < 0.0)) { 
	
	// If so, compute the frictional strength
	contact_law->computeFricStrength(tau_y[n], strength, n, it); 
	
	// Set 0 velocity in the normal direction
	for (UInt side = 0; side < 2; ++side) {
	  (*velocities[side])[n*dim+1] = 0.0;
	}
      }
      else{

	// Test if stresses are initially higher than initial strength
	
	// Direction X
	
	strength = nor_strength[n];
	
	if (tau_y[n] < strength){
	  (*velocities[0])[n*dim+1] = 0;
	  (*velocities[1])[n*dim+1] = 0;
	  (*intfc_trac)[n*dim+1] = tau_y[n];
	}

	else{
	  broken_nor = true; // Set the normal direction as being overshooted
	}
	
      }
      
      //velocities u_x & u_z (Mode II and III)
      for (UInt side = 0; side < 2; ++side) {

	strength = shr_strength[n];
	
	if (shear_stress[n] < strength) {
	  
	  for (UInt k = 0; k < 2; ++k) {
	    (*velocities[side])[n*dim+2*k] = 0.0;
	  }

	  // Set the tractions as the resistance
	  
	  (*intfc_trac)[n*dim+0] = tau_x[n];
	  (*intfc_trac)[n*dim+2] = tau_z[n];

	}
	
	else{
	  broken_shr = true; // Set the shear strength as being overshooted
	}
      }
      
      if ((broken_shr == true)||(broken_nor == true)){

	Real stress = total_stress[n];
	std::vector<Real> stress_ratio;
	if (stress > 0)
	  stress_ratio = {tau_x[n]/stress,tau_y[n]/stress,tau_z[n]/stress};
	else
	  stress_ratio = {0,0,0};	  

	Real radiation_damping = accoust / 2 * sqrt( (stress_ratio[0] * stress_ratio[0] + stress_ratio[2]*stress_ratio[2] + stress_ratio[1] * stress_ratio[1] * eta[0] * eta[0]) );
	
	Real rate = 0;
	
	// There is no opening at the first time step

	Real op_eq = 0;
	
	// If op_eq is larger than 1, the material is completely broken!
	  
	if ((op_eq>=1)||(nor_strength[n] * shr_strength[n] == 0)) {
	    
	  bool in_contact = ((nor_strength[n] == 0)&&(cRacklet::is_negative(0)));
	  // the case of contact is handled by the associated ContactLaw
	    
	  if (!in_contact) { 
	      
	    ind_crack[n] = 2;
	    nor_strength[n] = 0.0;
	    shr_strength[n] = 0.0;
	  }
	  
	  rate = stress / radiation_damping;

	  for (UInt k = 0; k < 3; ++k) {
	    (*intfc_trac)[n*dim+k] = 0;
	  }
	  
	}
	
	else{

	  // We need to compute the velocity for which we have equilibrium along the fault
	  
	  rate = (v_x_top[n] - v_x_bot[n]) * (v_x_top[n] - v_x_bot[n]) + (v_y_top[n] - v_y_bot[n]) * (v_y_top[n] - v_y_bot[n]) + (v_z_top[n] - v_z_bot[n])*(v_z_top[n] - v_z_bot[n]);
	  rate = sqrt(rate);
	  rate *= c_s;
	  
	  Real converg_v=1.0;
	  Real converg_s=1.0;
	  Real gamma = 0.1;
	  Real delta_v = 0.0;
	  Real delta_s,K;
	  Real m=0;
	  Real local_strength;
	  
	  while ((!cRacklet::has_converged(converg_v))&&(!cRacklet::has_converged(converg_s))) { //&&(is_good)) {
	    
	    rate += delta_v;
	    
	    // Compute the strength
	    
	    Real strength_v0 = max_nor_strength[n] * (1-op_eq);
	    
	    local_strength = formulation->getStrength(strength_v0,rate,lim_velocity[n]);
	    
	    Real tangent = formulation->getTangent(strength_v0,rate,lim_velocity[n]);
	    
	    // The radiation damping slope is an effective one that depends on the direction of the solicitation. For shear the slope is mu / 2 cs while for normal opening it is mu eta / 2 cs

	    delta_s = stress - local_strength - radiation_damping*rate;
	    
	    K = radiation_damping + tangent; 

	    delta_v = delta_s / K;
	    
	    /*if ((delta_v > 0))
	      delta_v = std::min(std::abs(delta_v),gamma*std::abs(rate));
	      else
	      delta_v = -std::min(std::abs(delta_v),gamma*std::abs(rate));
	    */	    
	    converg_v = delta_v/rate;
	    converg_s = delta_s;
	    m += 1;
	  }
	  
	  // Set to be in the cohesive zone and reduce the strength!
	  
	  ind_crack[n] = 1;
	  nor_strength[n] = max_nor_strength[n] * (1-op_eq);
	  shr_strength[n] = max_shr_strength[n] * (1-op_eq);
	}

	// Distribute the stresses and the velocities based on the stress ratio

	for (UInt s = 0; s < 3; ++s) {
	
	  (*velocities[0])[n*dim+s] = 0.5*(rate*stress_ratio[s])/c_s;
	  (*velocities[1])[n*dim+s] = -0.5*(rate*stress_ratio[s])/c_s;
	
	}
	(*intfc_trac)[n*dim+0] = stress*stress_ratio[0] - accoust*(rate*stress_ratio[0])/2;
	(*intfc_trac)[n*dim+2] = stress*stress_ratio[2] - accoust*(rate*stress_ratio[2])/2;
	(*intfc_trac)[n*dim+1] = stress*stress_ratio[1] - eta[0]*accoust*(rate*stress_ratio[1])/2;
	
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::updateInterfaceConditions() {

  Real strength;
  Real shr_trac;
  Real shr_velo;
  
  // Get the stresses
  
  CrackProfile tau_x = stresses[0]->getStridedPart(3,0);
  CrackProfile tau_y = stresses[0]->getStridedPart(3,1);
  CrackProfile tau_z = stresses[0]->getStridedPart(3,2);

  CrackProfile v_x_top = velocities[0]->getStridedPart(3,0);
  CrackProfile v_y_top = velocities[0]->getStridedPart(3,1);
  CrackProfile v_z_top = velocities[0]->getStridedPart(3,2);

  CrackProfile v_x_bot = velocities[1]->getStridedPart(3,0);
  CrackProfile v_y_bot = velocities[1]->getStridedPart(3,1);
  CrackProfile v_z_bot = velocities[1]->getStridedPart(3,2);
  
  CrackProfile * shr_opening = datas[_shear_displacement_jumps];
  CrackProfile * nor_opening = datas[_normal_displacement_jumps];
  
  CrackProfile shear_stress = tau_x*tau_x + tau_z*tau_z;
  shear_stress.squareRoot();

  CrackProfile total_stress = tau_y * tau_y + shear_stress * shear_stress;;
  total_stress.squareRoot();
  
  UInt n=0;
  for (UInt x = 0; x < n_ele[0]; ++x) {
    for (UInt z = 0; z < n_ele[1]; ++z) {
      n=x+z*n_ele[0];
      
      bool broken_nor = false;  
      bool broken_shr = false;  
      
      // Check if the element is broken AND in compression.
      if((this->nor_strength[n]==0)&&(tau_y[n] < 0.0)) { 
	
	// If so, compute the frictional strength
	contact_law->computeFricStrength(tau_y[n], strength, n, it); 
	
	// Set 0 velocity in the normal direction
	for (UInt side = 0; side < 2; ++side) {
	  (*velocities[side])[n*dim+1] = 0.0;
	}
      }
      else{

	// Test if stresses are initially higher than initial strength
	
	// Direction X
	
	strength = nor_strength[n];
	
	if (tau_y[n] <= strength){
	  (*velocities[0])[n*dim+1] = 0;
	  (*velocities[1])[n*dim+1] = 0;
	  (*intfc_trac)[n*dim+1] = tau_y[n];

	}

	else{
	  broken_nor = true; // Set the normal direction as being overshooted
	}
	
      }
      
      //velocities u_x & u_z (Mode II and III)
      for (UInt side = 0; side < 2; ++side) {

	strength = shr_strength[n];
	
	if (shear_stress[n] <= strength) {
	  
	  for (UInt k = 0; k < 2; ++k) {
	    (*velocities[side])[n*dim+2*k] = 0.0;
	  }

	  (*intfc_trac)[n*dim+0] = tau_x[n];
	  (*intfc_trac)[n*dim+2] = tau_z[n];

	}
	
	else{
	  broken_shr = true; // Set the shear strength as being overshooted
	}
      }
      
      if ((broken_shr == true)||(broken_nor == true)){

	Real stress = total_stress[n];
	std::vector<Real> stress_ratio;
	if (stress > 0)
	  stress_ratio = {tau_x[n]/stress,tau_y[n]/stress,tau_z[n]/stress};
	else
	  stress_ratio = {0,0,0};	  

	    Real radiation_damping = accoust / 2 * sqrt( (stress_ratio[0] * stress_ratio[0] + stress_ratio[2]*stress_ratio[2] + stress_ratio[1] * stress_ratio[1] * eta[0] * eta[0])) ;
	
	Real rate = 0;
	
	// Compute the equivalent opening
	    
	Real op_eq = sqrt(((*nor_opening)[n]/crit_nor_opening[n])*((*nor_opening)[n]/crit_nor_opening[n])+	   
			  ((*shr_opening)[n]/crit_shr_opening[n])*((*shr_opening)[n]/crit_shr_opening[n]));

	// If op_eq is larger than 1, the material is completely broken!
	  
	if ((op_eq>=1)||(nor_strength[n] * shr_strength[n] == 0)) {
	    
	  bool in_contact = ((nor_strength[n] == 0)&&(cRacklet::is_negative((*nor_opening)[n])));
	  // the case of contact is handled by the associated ContactLaw
	    
	  if (!in_contact) { 
	      
	    ind_crack[n] = 2;
	    nor_strength[n] = 0.0;
	    shr_strength[n] = 0.0;
	  }
	  
	  rate = stress / radiation_damping;

	  for (UInt k = 0; k < 3; ++k) {
	    (*intfc_trac)[n*dim+k] = 0;
	  }
	  
	}
	
	else{

	  // We need to compute the velocity for which we have equilibrium along the fault
	  
	  rate = (v_x_top[n] - v_x_bot[n]) * (v_x_top[n] - v_x_bot[n]) + (v_y_top[n] - v_y_bot[n]) * (v_y_top[n] - v_y_bot[n]) + (v_z_top[n] - v_z_bot[n])*(v_z_top[n] - v_z_bot[n]);
	  rate = sqrt(rate);
	  rate *= c_s;
	  
	  Real converg_v=1.0;
	  Real converg_s=1.0;
	  Real gamma = 0.1;
	  Real delta_v = 0.0;
	  Real delta_s,K;
	  Real m=0;
	  Real local_strength;

	  Real strength_v0 = max_nor_strength[n] * (1-op_eq);
	  
	  while ((!cRacklet::has_converged(converg_v))&&(!cRacklet::has_converged(converg_s))) { //&&(is_good)) {
	  
	    rate += delta_v;
	    
	    // Exponent that depends on the particle velocity (rate)

	    local_strength = formulation->getStrength(strength_v0,rate,lim_velocity[n]);
	    
	    Real tangent = formulation->getTangent(strength_v0,rate,lim_velocity[n]);
	    
	    K = radiation_damping + tangent; 

	    delta_s = stress - local_strength - radiation_damping*rate;
	    
	    delta_v = delta_s / K;
	    
	    if ((delta_v > 0))
	      delta_v = std::min(std::abs(delta_v),gamma*std::abs(lim_velocity[n]));
	    /*
	      else
	      delta_v = -std::min(std::abs(delta_v),gamma*std::abs(lim_velocity[n]));
	    */
	    converg_v = delta_v/rate;
	    converg_s = delta_s/stress;
	    m += 1;

	    if(m > 100){
	      std::cout<< "Iterations " << m << std::endl;
	      std::cout<< "Tangent " << tangent << std::endl;
	      std::cout<< " K" << K << std::endl;
	      std::cout<< "delta_v : " << delta_v << std::endl;
	      std::cout<< "delta_s :" << delta_s << std::endl;
	      std::cout<< "converg_v:  " << converg_v << std::endl;
	      std::cout<< "converg_s: " << converg_s << std::endl;
	    
	    }

	    if(m>20){
	      cRacklet::error("Solution has not been found");
	      break;
	    }
	  }
	  
	  // Set to be in the cohesive zone and reduce the strength!
	  
	  ind_crack[n] = 1;
	  nor_strength[n] = max_nor_strength[n] * (1-op_eq);
	  shr_strength[n] = max_shr_strength[n] * (1-op_eq);
	}

	// Distribute the stresses and the velocities based on the stress ratio

	for (UInt s = 0; s < 3; ++s) {
	
	  (*velocities[0])[n*dim+s] = 0.5*(rate*stress_ratio[s])/c_s;
	  (*velocities[1])[n*dim+s] = -0.5*(rate*stress_ratio[s])/c_s;
	}

	(*intfc_trac)[n*dim+0] = stress*stress_ratio[0] - accoust*(rate*stress_ratio[0])/2;
	(*intfc_trac)[n*dim+2] = stress*stress_ratio[2] - accoust*(rate*stress_ratio[2])/2;
	(*intfc_trac)[n*dim+1] = stress*stress_ratio[1] - eta[0]*accoust*(rate*stress_ratio[1])/2;

      }
    }
  }
}
