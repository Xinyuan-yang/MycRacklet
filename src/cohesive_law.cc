/* -------------------------------------------------------------------------- */
#include "cohesive_law.hh"
#include <fstream>
#include <algorithm>
/* -------------------------------------------------------------------------- */
void CohesiveLaw::initInterfaceConditions() {
  
  memcpy(&max_nor_strength[0],&nor_strength[0],nor_strength.size()*sizeof(Real));
  memcpy(&max_shr_strength[0],&shr_strength[0],shr_strength.size()*sizeof(Real));
  
  computeInitialVelocities();
}

/* -------------------------------------------------------------------------- */
void CohesiveLaw::preventSurfaceOverlapping(std::shared_ptr<ContactLaw> contactlaw) {

  allow_overlapping = false;
  contact_law = contactlaw;  
}

/* -------------------------------------------------------------------------- */
void CohesiveLaw::restart(bool pausing, UInt nele_2d) {
  
  if(contact_law)
    contact_law->restart(pausing,nele_2d);
}

/* -------------------------------------------------------------------------- */
void CohesiveLaw::updateInterfaceConditions() {
  updateCohesiveLaw();
  computeVelocities();
}

/* -------------------------------------------------------------------------- */
void CohesiveLaw::computeInitialVelocities() {

  Real strength;
  Real shr_trac;
  Real shr_velo;
  std::vector<Real> temp_f(2);

  std::vector<CrackProfile*> loads = {datas[_top_loading],datas[_bottom_loading]};
  
  UInt i=0;
  for (UInt h = 0; h < n_ele[0]; ++h) {
    for (UInt j = 0; j < n_ele[1]; ++j) {
      i=h+j*n_ele[0];
    
      if((nor_strength[i]==0)&&((*loads[0])[i*dim+1] < 0.0)) { 
      
	contact_law->computeFricStrength((*loads[0])[i*dim+1], strength, i, it); 
      
	for (UInt side = 0; side < 2; ++side) {
	  (*velocities[side])[i*dim+1] = 0.0;
	}
      }
      else{//velocities u2
	strength = shr_strength[i];
	(*velocities[0])[i*dim+1] = std::max(((*loads[0])[i*dim+1]-nor_strength[i])/(mu[0]*eta[0]),0.0);
	(*velocities[1])[i*dim+1] = std::min((zeta/ksi)*(nor_strength[i]-(*loads[1])[i*dim+1])/(mu[0]*eta[1]),0.0);
      }

      //velocities u1 & u3
      for (UInt side = 0; side < 2; ++side) {
  
	for (UInt k = 0; k < 2; ++k) {
	  temp_f[k] = (*loads[side])[i*dim+2*k];
	}
      
	shr_trac = sqrt(temp_f[0]*temp_f[0]+temp_f[1]*temp_f[1]);
	if (shr_trac ==0) {
	
	  for (UInt k = 0; k < 2; ++k) {
	  
	    (*velocities[side])[i*dim+2*k] = 0.0;
	  }
	}
	else {
	  if (side==0) shr_velo = std::max((shr_trac-strength)/mu[0],0.0);
	  else shr_velo = std::min((zeta/ksi)*(strength - shr_trac)/mu[0],0.0);
	  
	  for (UInt k = 0; k < 2; ++k) {
	    (*velocities[side])[i*dim+2*k] = shr_velo*temp_f[k]/shr_trac;
	  } 
	}
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveLaw::updateCohesiveLaw() {

  CrackProfile * shr_opening = datas[_shear_displacement_jumps];
  CrackProfile * nor_opening = datas[_normal_displacement_jumps];

  UInt n_ele = max_nor_strength.size();
  Real aux;

  for (UInt i = 0; i < n_ele; ++i) {     

    if (((*nor_opening)[i]==0)&&((*shr_opening)[i]==0)&&(ind_crack[i]!=2)){

      nor_strength[i] = max_nor_strength[i];
      shr_strength[i] = max_shr_strength[i];

    }
    
    else {
      
      aux = sqrt(((*nor_opening)[i]/crit_nor_opening[i])*((*nor_opening)[i]/crit_nor_opening[i])+
		 ((*shr_opening)[i]/crit_shr_opening[i])*((*shr_opening)[i]/crit_shr_opening[i]));
      
      if ((aux>=1)||(nor_strength[i] * shr_strength[i] == 0)) {

	bool in_contact = ((nor_strength[i] == 0)&&(cRacklet::is_negative((*nor_opening)[i])));
	// the case of contact is handled by the associated ContactLaw
	
	if (!in_contact) { 

	  ind_crack[i] = 2;
	  nor_strength[i] = 0.0;
	  shr_strength[i] = 0.0;
	}
      }

      else {

	ind_crack[i] = 1;
	nor_strength[i] = max_nor_strength[i] * (1-aux);
	shr_strength[i] = max_shr_strength[i] * (1-aux);
      }
    }
  }    
}

/* -------------------------------------------------------------------------- */
void CohesiveLaw::computeVelocities(){

  CrackProfile deltaStresses(n_ele,dim);
  std::vector<Real> temp_veloc(dim);
  Real trac;
   
  deltaStresses = (*stresses[0]) - (*stresses[1]);
   
  Real cste = 1/(mu[0]*(1+ksi/zeta)); 

  (*velocities[0]) =  deltaStresses * cste;

  for (UInt i = 0; i < n_ele[0]; ++i) {
    for (UInt j = 0; j < n_ele[1]; ++j) {
      (*velocities[0])[(i*dim+1)+j*n_ele[0]*dim] *= (1+ksi/zeta)/(eta[0]+ksi*eta[1]/zeta); 
    }
  }
  (*velocities[1])=(*velocities[0]); 

  for (UInt i = 0; i < n_ele[0]; ++i) {
    for (UInt j = 0; j < n_ele[1]; ++j) {

      trac = (*stresses[0])[(i*dim+1)+j*n_ele[0]*dim] - mu[0]*eta[0]* (*velocities[0])[(i*dim+1)+j*n_ele[0]*dim];
      if ((nor_strength[i+n_ele[0]*j] < trac)||(nor_strength[i+n_ele[0]*j]==0)) computeIndepNormalVelocities(i,j);
      else {
	(*intfc_trac)[(i*dim+1)+j*n_ele[0]*dim] = trac;
	computeShearVelocities(shr_strength[i+n_ele[0]*j], i+j*n_ele[0]);
      }
    }
  }
}

