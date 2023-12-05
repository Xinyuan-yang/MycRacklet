/**
 * @file   cohesive_law.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Sun Jan  6 19:35:47 2013
 *
 * @brief  Implementation of the CohesiveLaw class
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

  std::vector<Real> * nor_strength = datas[_normal_strength];
  DataRegister::restartData(*nor_strength,"restart_normal_strength.cra",pausing, nele_2d);
  std::vector<Real> * shr_strength = datas[_shear_strength];
  DataRegister::restartData(*shr_strength,"restart_shear_strength.cra",pausing, nele_2d);
  
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
      
      if((nor_strength[i]==res_nor_strength[i])&&((*loads[0])[i*dim+1] < 0.0)) { 
	
	contact_law->computeFricStrength((*loads[0])[i*dim+1], strength, i, it); 
	
	for (UInt side = 0; side < 2; ++side) {
	  (*velocities[side])[i*dim+1] = 0.0;
	}
      }
      else{//velocities u2
	(*velocities[0])[i*dim+1] = std::max(((*loads[0])[i*dim+1]-nor_strength[i])/(mu[0]*eta[0])*cs[0],0.0);
	(*velocities[1])[i*dim+1] = std::min((nor_strength[i]-(*loads[1])[i*dim+1])/(mu[1]*eta[1])*cs[1],0.0);
      }
      
      //velocities u1 & u3
      for (UInt side = 0; side < 2; ++side) {
	strength = shr_strength[i];
	
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
	  if (side==0) shr_velo = std::max((shr_trac-strength)/mu[0]*cs[0],0.0);
	  else shr_velo = std::min((strength - shr_trac)/mu[0]*cs[1],0.0);
	  
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
      
      if ((aux>=1)||(nor_strength[i]==res_nor_strength[i])||(shr_strength[i]==res_shr_strength[i])) {

	bool in_contact = ((nor_strength[i] == res_nor_strength[i])&&(cRacklet::is_negative((*nor_opening)[i])));
	// the case of contact is handled by the associated ContactLaw
	
	if (!in_contact) { 

	  ind_crack[i] = 2;
	  nor_strength[i] = res_nor_strength[i];
	  shr_strength[i] = res_shr_strength[i];
	}
      }

      else {
	ind_crack[i] = 1;
	nor_strength[i] = max_nor_strength[i] - (max_nor_strength[i]-res_nor_strength[i])*aux;
	shr_strength[i] = max_shr_strength[i] - (max_shr_strength[i]-res_shr_strength[i])*aux;

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
   
  Real cste = 1/(mu[0]/cs[0]+mu[1]/cs[1]); 
  
  (*velocities[0]) =  deltaStresses * cste;

  // Correction for the normal velocity
  for (UInt i = 0; i < n_ele[0]; ++i) {
    for (UInt j = 0; j < n_ele[1]; ++j) {
      (*velocities[0])[(i*dim+1)+j*n_ele[0]*dim] *= 1/(mu[0]*eta[0]/cs[0]+mu[1]*eta[1]/cs[1])*(1/cste); 
    }
  }
  (*velocities[1])=(*velocities[0]); 
  
  for (UInt i = 0; i < n_ele[0]; ++i) {
    for (UInt j = 0; j < n_ele[1]; ++j) {
      
      trac = (*stresses[0])[(i*dim+1)+j*n_ele[0]*dim] - mu[0]*eta[0]* (*velocities[0])[(i*dim+1)+j*n_ele[0]*dim] / cs[0];
      
      if ((nor_strength[i+n_ele[0]*j] < trac)||(nor_strength[i+n_ele[0]*j]==0))
	computeIndepNormalVelocities(i,j);
      else {
	(*intfc_trac)[(i*dim+1)+j*n_ele[0]*dim] = trac;
	computeShearVelocities(shr_strength[i+n_ele[0]*j], i+j*n_ele[0]);
      }
    }
  }
}


/* -------------------------------------------------------------------------- */
void CohesiveLaw::correctVelocities(std::vector<Real> vel_correction){
  
  for (UInt i = 0; i < n_ele[0]; ++i) {
    for (UInt j = 0; j < n_ele[1]; ++j) {
      for (UInt d = 0; d < dim; ++d) {
	(*velocities[0])[(i+j*n_ele[0])*dim+d] += vel_correction[(i+j*n_ele[0])*dim+d];
      }
    }
  }

  // Change also the bottom velocities
  (*velocities[1])=(*velocities[0]); 
  
}
/* -------------------------------------------------------------------------- */
