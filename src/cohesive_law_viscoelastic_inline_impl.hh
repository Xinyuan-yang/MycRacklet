/**
 * @file   cohesive_law_viscoelastic_inline.cc
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Wed Nov 18 18:56:06 2020
 *
 * @brief  Implementation of the inline functions of the CohesiveLawViscoelastic class
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
void CohesiveLawViscoelastic::computeIndepNormalVelocities(UInt ix, UInt iz){
  
  std::vector<Real> cmpted_stress(2);
  Real delta_overlap;

  UInt i = ix+iz*n_ele[0];

  for (UInt side = 0; side < 2; ++side) {
    cmpted_stress[side] = (*stresses[side])[i*dim+1]; 
  }

  // Compute the velocity for the top part. Here we assume it's the same material!
    
  Real rate = this->prev_nor_vel[i];
  
  Real converg_v=1.0;
  Real converg_s=1.0;
  Real gamma = 0.1;
  Real delta_v = 0.0;
  Real delta_s,K;
  Real m=0;
  Real local_strength,tangent;
  
  while ((!cRacklet::has_converged(converg_v))&&(!cRacklet::has_converged(converg_s))) { //&&(is_good)) {
    
    rate += delta_v;
    
    // Compute the strength
    
    local_strength = formulation->getStrength(max_nor_strength[i] * (1-this->op_eq[i]),rate,lim_velocity[i]);
    
    tangent = formulation->getTangent(max_nor_strength[i] * (1-this->op_eq[i]),rate,lim_velocity[i]);
    
    // The radiation damping slope is mu eta / 2 cs
    
    delta_s = cmpted_stress[0]- local_strength - mu[0]*eta[0]*rate/c_s/2;
    
    K = mu[0]*eta[0]/c_s/2 + tangent; 
    
    delta_v = delta_s / K;
    
    /*if ((delta_v > 0))
      delta_v = std::min(std::abs(delta_v),gamma*std::abs(lim_velocity[i]));
    else
      delta_v = -std::min(std::abs(delta_v),gamma*std::abs(lim_velocity[i]));
    */
    converg_v = delta_v/rate;
    converg_s = delta_s/std::max(1.,cmpted_stress[0]);
    m += 1;

    if(m>20){
      std::cout << "Normal..." << std::endl;
      std::cout << "Delta v : " << delta_v << std::endl;
      std::cout << "Delta s : " << delta_s << std::endl;
      std::cout << "Error v : " << converg_v << std::endl;
      std::cout << "Error s : " << converg_s << std::endl;
      cRacklet::error("Solution has not been found");
    }
    
  }

  Real dt = dxmin*beta/(std::max(cs[0],cs[1]));
  
  delta_overlap = (*displacements[0])[i*dim+1] - (*displacements[1])[i*dim+1] + dt*rate; 
  
  if (cRacklet::is_negative(delta_overlap)&&(!allow_overlapping)){
    computeContactVelocities(ix, iz);
  }
  else {

    this->prev_nor_vel[i] = rate;
    
    (*velocities[0])[i*dim+1] = 0.5*rate;
    (*velocities[1])[i*dim+1] = -0.5*rate;
    (*intfc_trac)[i*dim+1] = cmpted_stress[0]-mu[0]*eta[0]*rate/c_s/2;       
    computeShearVelocities(shr_strength[i], i);
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::computeContactVelocities(UInt ix, UInt iz){
  
  Real aux;
  Real temp_velot;
  Real temp_trac;
  Real strength;
  std::vector<Real> cmpted_stress(2);
  
  UInt i = ix+iz*n_ele[0];
  
  for (UInt side = 0; side < 2; ++side) {
    cmpted_stress[side] = (*stresses[side])[i*dim+1];
  }
  
  aux = ((*displacements[0])[i*dim+1] - (*displacements[1])[i*dim+1])/(dxmin*beta);
  temp_velot = 1/(eta[0]/cs[0]+eta[1]/cs[1]/(mu[0]/mu[1]))*((cmpted_stress[0]-cmpted_stress[1])/mu[0] - aux*eta[1]/cs[1]/(mu[0]/mu[1]));
  (*velocities[0])[i*dim+1] = temp_velot;
  (*velocities[1])[i*dim+1] = temp_velot + aux;
  
  temp_trac = cmpted_stress[0] - eta[0]*mu[0]*temp_velot/cs[0];
  (*intfc_trac)[i*dim+1] = temp_trac;
  
  contact_law->computeFricStrength(temp_trac, strength, i, it);
  
  computeShearVelocities(strength, i);
  
  ind_crack[i] = 3;
  fric_strength[i] = strength;
}

/* -------------------------------------------------------------------------- */
void CohesiveLawViscoelastic::computeShearVelocities(Real strength, UInt i) {
  
  std::vector<Real> cmpted_stress(2);
  Real shr_trac;
  
  // Look only at the top side...
  
  for (UInt j = 0; j < 2; ++j) {
    cmpted_stress[j] = (*stresses[0])[i*dim+2*j];
  }
  
  Real rate = this->prev_shr_vel[i];
  
  shr_trac = sqrt((cmpted_stress[0]*cmpted_stress[0])+(cmpted_stress[1]*cmpted_stress[1])); 
  
  Real local_strength = formulation->getStrength(max_shr_strength[i] * (1-this->op_eq[i]),rate,lim_velocity[i]);
  
  if (shr_trac <= local_strength){
    for (UInt j = 0; j < (dim-1); ++j) {
      (*intfc_trac)[i*dim+2*j] = cmpted_stress[j];
    }
  }else if(shr_trac == 0){
    for (UInt j = 0; j < (dim-1); ++j) {
      (*intfc_trac)[i*dim+2*j] = cmpted_stress[j];
      (*velocities[0])[i*dim+2*j]=0;
      (*velocities[1])[i*dim+2*j]=0;
    }
  }
  else{
    
    Real converg_v=1.0;
    Real converg_s=1.0;
    Real gamma = 0.1;
    Real delta_v = 0.0;
    Real delta_s,K;
    Real m=0;
    Real tangent;
    
    while ((!cRacklet::has_converged(converg_v))&&(!cRacklet::has_converged(converg_s))) { //&&(is_good)) {
    
      rate += delta_v;
    
      // Compute the strength
      
      local_strength = formulation->getStrength(max_shr_strength[i] * (1-op_eq[i]),rate,lim_velocity[i]);
      
      tangent = formulation->getTangent(max_shr_strength[i] * (1-op_eq[i]),rate,lim_velocity[i]);
      
      // The radiation damping slope is mu / 2 cs
      
      delta_s = shr_trac - local_strength - mu[0]*rate/c_s/2;
    
      K = mu[0]/c_s/2 + tangent; 
      
      delta_v = delta_s / K;
      
      /*if ((delta_v > 0))
	delta_v = std::min(std::abs(delta_v),gamma*std::abs(lim_velocity[i]));
	else
	delta_v = -std::min(std::abs(delta_v),gamma*std::abs(lim_velocity[i]));
      */
      converg_v = delta_v/rate;
      converg_s = delta_s/shr_trac;
      m += 1;
      
      if(m>20){
	std::cout << "Delta v : " << delta_v << std::endl;
	std::cout << "Delta s : " << delta_s << std::endl;
	std::cout << "Error v : " << converg_v << std::endl;
	std::cout << "Error s : " << converg_s << std::endl;
	cRacklet::error("Solution has not been found");
      }
      
    }
    
    for (UInt j = 0; j < 2; ++j) {
      (*velocities[0])[i*dim+2*j] = 0.5*rate*cmpted_stress[j]/shr_trac;
      (*velocities[1])[i*dim+2*j] = -0.5*rate*cmpted_stress[j]/shr_trac;    
      (*intfc_trac)[i*dim+2*j] = local_strength*cmpted_stress[j]/shr_trac;
      
    }

    this->prev_shr_vel[i] = rate;

  }
}
