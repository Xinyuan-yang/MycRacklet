/**
 * @file   cohesive_law_viscoelastic.cc
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Wed Nov 18 18:56:06 2020
 *
 * @brief  Implementation of the CohesiveLawViscoelastic class
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
#include "cohesive_law_viscoelastic.hh"
#include <fstream>
#include <algorithm>
#include <math.h> 
/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::initLinearFormulation() {
  formulation = std::make_shared<ViscoelasticFormulation>();
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::initQuadraticFormulation() {
  formulation = std::make_shared<ViscoelasticQuadraticFormulation>();
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::initPowerLawFormulation() {
  formulation = std::make_shared<ViscoelasticPowerLawFormulation>();
  
  memcpy(&this->lim_velocity[0],&this->lim_velocity[0],this->lim_velocity.size()*sizeof(Real));
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::initInterfaceConditions() {
  
  memcpy(&this->max_nor_strength[0],&this->nor_strength[0],this->nor_strength.size()*sizeof(Real));
  memcpy(&this->max_shr_strength[0],&this->shr_strength[0],this->shr_strength.size()*sizeof(Real));
  
  this->computeInitialVelocities();
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::updateInterfaceConditions() {
  this->updateCohesiveLaw();
  this->computeVelocities();
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
  
  // Get the stresses

  std::vector<CrackProfile*> loads = {datas[_top_loading],datas[_bottom_loading]};

  for (UInt side = 0; side < 2; ++side) {
    for (UInt i = 0; i < n_ele[0]; ++i) {
      for (UInt j = 0; j < n_ele[1]; ++j) {
	
	(*stresses[side])[(i*dim+1)+j*n_ele[0]*dim] = (*loads[side])[(i*dim+1)+j*n_ele[0]*dim];
      }
    }
  }

  this->computeVelocities();
  
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::updateCohesiveLaw() {

  CrackProfile * shr_opening = datas[_shear_displacement_jumps];
  CrackProfile * nor_opening = datas[_normal_displacement_jumps];
  
  UInt n_ele = max_nor_strength.size();
  
  for (UInt i = 0; i < n_ele; ++i) {     
    
    if (((*nor_opening)[i]==0)&&((*shr_opening)[i]==0)&&(ind_crack[i]!=2)){
      
      nor_strength[i] = max_nor_strength[i];
      shr_strength[i] = max_shr_strength[i];
      
    }
    
    else {
      
      this->op_eq[i] = std::min(sqrt(((*nor_opening)[i]/crit_nor_opening[i])*((*nor_opening)[i]/crit_nor_opening[i])+
				     ((*shr_opening)[i]/crit_shr_opening[i])*((*shr_opening)[i]/crit_shr_opening[i])),1.);
      
      
      
      if ((this->op_eq[i]>=1)||(nor_strength[i] * shr_strength[i] == 0)) {

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
	
	Real rate_nor = ((*velocities[0])[i*dim+1]-(*velocities[1])[i*dim+1]);
	Real rate_shr = sqrt( ( (*velocities[0])[i*dim+0]-(*velocities[1])[i*dim+0])*( (*velocities[0])[i*dim+0]-(*velocities[1])[i*dim+0]) + ( (*velocities[0])[i*dim+2]-(*velocities[1])[i*dim+2]) * ( (*velocities[0])[i*dim+2]-(*velocities[1])[i*dim+2]));
	
	nor_strength[i] = formulation->getStrength(max_nor_strength[i] * (1-this->op_eq[i]),rate_nor,lim_velocity[i]);
	shr_strength[i] = formulation->getStrength(max_shr_strength[i] * (1-this->op_eq[i]),rate_shr,lim_velocity[i]);;
      }
    }
  }    
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::computeVelocities(){

  CrackProfile deltaStresses(n_ele,dim);
  std::vector<Real> temp_veloc(dim);
  Real trac;
   
  deltaStresses = (*stresses[0]) - (*stresses[1]);
   
  Real cste = 1/(mu[0]/cs[0]+mu[1]/cs[1]); 

  (*velocities[0]) =  deltaStresses * cste;

  for (UInt i = 0; i < n_ele[0]; ++i) {
    for (UInt j = 0; j < n_ele[1]; ++j) {
      (*velocities[0])[(i*dim+1)+j*n_ele[0]*dim] *= 1/(mu[0]*eta[0]/cs[0]+mu[1]*eta[1]/cs[1])*(1/cste); 
    }
  }
  (*velocities[1])=(*velocities[0]); 
  
  for (UInt i = 0; i < n_ele[0]; ++i) {
    for (UInt j = 0; j < n_ele[1]; ++j) {
      
      Real trac = (*stresses[0])[(i*dim+1)+j*n_ele[0]*dim] - mu[0]*eta[0]* (*velocities[0])[(i*dim+1)+j*n_ele[0]*dim]/cs[0];
      if ((nor_strength[i+n_ele[0]*j] < trac)||(nor_strength[i+n_ele[0]*j]==0))
	computeIndepNormalVelocities(i,j);
      else {
	(*intfc_trac)[(i*dim+1)+j*n_ele[0]*dim] = trac;
	computeShearVelocities(shr_strength[i+n_ele[0]*j], i+j*n_ele[0]);
      }
    }
  }
}

