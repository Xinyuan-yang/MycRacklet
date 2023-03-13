/**
 * @file   cohesive_law_inline_impl.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Sun Jan  6 19:35:47 2013
 *
 * @brief  Implementation of the inline functions of the CohesiveLaw class
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
void CohesiveLaw::computeIndepNormalVelocities(UInt ix, UInt iz){

  std::vector<Real> temp_veloc(2);
  std::vector<Real> cmpted_stress(2);
  Real delta_overlap;

  UInt i = ix+iz*n_ele[0];

  for (UInt side = 0; side < 2; ++side) {
    cmpted_stress[side] = (*stresses[side])[i*dim+1]; 
  }
  
  temp_veloc[0] = cs[0]/(mu[0]*eta[0]) * (cmpted_stress[0]-nor_strength[i]); 
  temp_veloc[1] = cs[1]/(mu[1]*eta[1]) * (nor_strength[i] - cmpted_stress[1]);

  Real dt = dxmin*beta/(std::max(cs[0],cs[1]));
  
  delta_overlap = (*displacements[0])[i*dim+1] - (*displacements[1])[i*dim+1] + dt*(temp_veloc[0]-temp_veloc[1]); 

  if (cRacklet::is_negative(delta_overlap)&&(!allow_overlapping)) computeContactVelocities(ix, iz); 
  else {

    for (UInt side = 0; side < 2; ++side) {

      (*velocities[side])[i*dim+1] = temp_veloc[side];
    }
    (*intfc_trac)[i*dim+1] = nor_strength[i]; 
    computeShearVelocities(shr_strength[i], i);
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveLaw::computeContactVelocities(UInt ix, UInt iz){
  
  Real aux;
  Real temp_velot;
  Real temp_trac;
  Real strength;
  std::vector<Real> cmpted_stress(2);

  Real dt = dxmin*beta/(std::max(cs[0],cs[1]));
  
  UInt i = ix+iz*n_ele[0];
  
  for (UInt side = 0; side < 2; ++side) {
    cmpted_stress[side] = (*stresses[side])[i*dim+1];
  }
  
  aux = ((*displacements[0])[i*dim+1] - (*displacements[1])[i*dim+1])/dt;
  temp_velot = 1/(eta[0]/cs[0]+eta[1]/cs[1]/(mu[0]/mu[1]))*((cmpted_stress[0]-cmpted_stress[1])/mu[0] - aux*eta[1]/cs[1]/(mu[0]/mu[1]));
  (*velocities[0])[i*dim+1] = temp_velot;
  (*velocities[1])[i*dim+1] = temp_velot + aux;
  
  temp_trac = cmpted_stress[0] - eta[0]*mu[0]*temp_velot/cs[0];
  (*intfc_trac)[i*dim+1] = temp_trac;
  
  if(contact_law!=NULL){
    contact_law->computeFricStrength(temp_trac, strength, i, it);
  }

  computeShearVelocities(strength, i);
  
  ind_crack[i] = 3;
  fric_strength[i] = strength;
}

/* -------------------------------------------------------------------------- */
void CohesiveLaw::computeShearVelocities(Real strength, UInt i) {
  
  std::vector<Real> trac(2);
  Real shr_trac;
  
  for (UInt j = 0; j < 2; ++j) {
    
    trac[j] = (*stresses[0])[i*dim+2*j] - mu[0]*(*velocities[0])[i*dim+2*j]/cs[0]; 
  }
  
  shr_trac = sqrt((trac[0]*trac[0])+(trac[1]*trac[1])); 

  if ((strength < shr_trac)||(strength==0)) computeIndepShearVelocities(strength, i);
  else{
    
    for (UInt j = 0; j < (dim-1); ++j) {
      
      (*intfc_trac)[i*dim+2*j] = trac[j];
    }
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveLaw::computeIndepShearVelocities(Real strength, UInt i){

  std::vector<Real> cmpted_stress(2);
  Real dyn_stress;
  Real shr_veloc;

  for (UInt side = 0; side < 2; ++side) {
     
    for (UInt j = 0; j < 2; ++j) {
       
      cmpted_stress[j] = (*stresses[side])[i*dim+2*j];
    }
    
    dyn_stress = sqrt((cmpted_stress[0]*cmpted_stress[0])+(cmpted_stress[1]*cmpted_stress[1])); 
    
    if (side==0)
      shr_veloc = cs[0]/mu[0]*(dyn_stress-strength); 
    
    else
      shr_veloc = cs[1]/mu[1]*(strength-dyn_stress);   

    for (UInt j = 0; j < 2; ++j) {
       
      if(dyn_stress==0)
	{(*velocities[side])[i*dim+2*j]=0;}
      else
	{(*velocities[side])[i*dim+2*j] = shr_veloc*cmpted_stress[j]/dyn_stress;}
      if (side==0){
	if(dyn_stress==0)
	  {(*intfc_trac)[i*dim+2*j] =0;}
	else
	  {(*intfc_trac)[i*dim+2*j] = strength*cmpted_stress[j]/dyn_stress;}
      }
    }    
  }
}
