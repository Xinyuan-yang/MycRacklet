/**
 * @file   cohesive_law_viscoelastic.hh
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Wed Nov 18 18:56:06 2020
 *
 * @brief  Class describing a coupled slip/opening weakening interface with velocity dependance
 *
 * @section LICENSE
 *
 * cRacklet - A spectral boundary integral method for interface fracture simulation
 * Copyright (©) 2012 - 2013 Fabian Barras
 *               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef __COHESIVE_LAW_VISCOELASTIC__
#define __COHESIVE_LAW_VISCOELASTIC__
/* -------------------------------------------------------------------------- */
#include "cohesive_law.hh"
/* -------------------------------------------------------------------------- */

class CohesiveLawViscoelastic : public CohesiveLaw {

  #include "cohesive_law_viscoelastic_formulations.hh"

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  CohesiveLawViscoelastic() 
    :CohesiveLaw() {
    
    UInt total_n_ele = n_ele[0]*n_ele[1];
    
    lim_velocity.resize(total_n_ele);
    op_eq.resize(total_n_ele);
    prev_nor_vel.resize(total_n_ele);
    prev_shr_vel.resize(total_n_ele);
    
    this->registerData(_lim_velocity, &lim_velocity);
    
    Real mu_top = this->getParameter<Real>("shear modulus top");
    Real mu_bot = this->getParameter<Real>("shear modulus bottom");
    Real cs_top = this->getParameter<Real>("shear wave speed top");
    Real cs_bot = this->getParameter<Real>("shear wave speed bottom");
    
    if ((mu_top==mu_bot)&&(cs_top==cs_bot)) {
      c_s = cs_top;
    }
    else
      cRacklet::error("Viscoelastic is only implemented for homogeneous properties");

  };

  virtual ~CohesiveLawViscoelastic();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /** formulation with linear strenghtening */
  void initLinearFormulation();
  /** formulation with power law strenghtening */
  void initPowerLawFormulation();
  
  /** Initialize interface fields */
  void initInterfaceConditions();
  /** Compute interface tractions and velocities in function of the new strength profile */
  void updateInterfaceConditions();
  /** Method used in restart framework but no history-dependant variable within this law */
  void restart(bool pausing=false, UInt nele_2d=0);
  
protected:
  
  /** compute velocities at t=0 */
  void computeInitialVelocities();
  /** update interface strength from cohesive law */
  void updateCohesiveLaw();
  /** Compute velocities */
  void computeVelocities();
  /** Compute normal velocities in case of relative slip */
  inline void computeIndepNormalVelocities(UInt ix, UInt iz);
  /** Compute shear velocities with a given shear strength */
  inline void computeShearVelocities(Real strength, UInt elem);
  /** Compute velocities in the case of contact at crack type */
  inline void computeContactVelocities(UInt ix, UInt iz);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /** Parameters */
  Real c_s;
  
  /** Limiting velocity */
  std::vector<Real> lim_velocity;

  /** Previous velocity used in computing the following one! */
  std::vector<Real> prev_nor_vel;
  std::vector<Real> prev_shr_vel;
  
  /** Equivalent opening */
  std::vector<Real> op_eq;
  
  // Abstract object representing the associated viscoelastic formulation
  std::shared_ptr<ViscoelasticFormulation> formulation;
  
};


/* -------------------------------------------------------------------------- */
/* inline functions                                       

                    */
/* -------------------------------------------------------------------------- */

inline CohesiveLawViscoelastic::~CohesiveLawViscoelastic(){
  
}							


#include "cohesive_law_viscoelastic_inline_impl.cc"

#endif /* __COHESIVE_LAW_VISCOELASTIC__ */
