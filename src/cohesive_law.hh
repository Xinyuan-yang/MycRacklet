/**
 * @file   cohesive_law.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Sun Jan  6 19:35:47 2013
 *
 * @brief  Class describing a coupled slip/opening weakening cohesive interface
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
#ifndef __COHESIVE_LAW__
#define __COHESIVE_LAW__
/* -------------------------------------------------------------------------- */
#include "cRacklet_common.hh"
#include "interface_law.hh"
#include "contact_law.hh"
#include <vector>
#include <math.h>
/* -------------------------------------------------------------------------- */

class CohesiveLaw : public InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  CohesiveLaw() 
    :InterfaceLaw() {

    UInt total_n_ele = n_ele[0]*n_ele[1];
    
    nor_strength.resize(total_n_ele);
    shr_strength.resize(total_n_ele);
    fric_strength.resize(total_n_ele);
    ind_crack.resize(total_n_ele);

    max_nor_strength.resize(total_n_ele);
    max_shr_strength.resize(total_n_ele);
    crit_nor_opening.resize(total_n_ele);
    crit_shr_opening.resize(total_n_ele);
    
    this->registerData(_normal_strength, &nor_strength);
    this->registerData(_shear_strength, &shr_strength);
    this->registerData(_frictional_strength, &fric_strength);
    this->registerData(_id_crack, &ind_crack);
    this->registerData(_critical_normal_opening, &crit_nor_opening);
    this->registerData(_critical_shear_opening, &crit_shr_opening); 
    
    mu = {getParameter<Real>("shear modulus top"),getParameter<Real>("shear modulus bottom")};
    stresses = {datas[_top_dynamic_stress],datas[_bottom_dynamic_stress]};
    velocities = {datas[_top_velocities],datas[_bottom_velocities]};
    displacements = {datas[_top_displacements],datas[_bottom_displacements]};
    intfc_trac = datas[_interface_tractions];

    contact_law = NULL;
    allow_overlapping = true;
  };

  virtual ~CohesiveLaw();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  
  /** Initialize interface fields */
  void initInterfaceConditions();
  /** Compute interface tractions and velocities in function of the new strength profile */
  void updateInterfaceConditions();
  /** Method used in restart framework but no history-dependant variable within this law */
  void restart(bool pausing=false, UInt nele_2d=0);
  // Prevent the overlapping of the two surfaces and associate a ContactLaw
  void preventSurfaceOverlapping(std::shared_ptr<ContactLaw> contactlaw);

  // Count the number of elements being part of the process zone (ind_crack = 1)
  UInt getNbCohesiveNodes(){return std::count(this->ind_crack.begin(),this->ind_crack.end(),1);};
  // Count the number of elements being part of the process zone (ind_crack = 2)
  UInt getNbBrokenNodes(){return std::count(this->ind_crack.begin(),this->ind_crack.end(),2);};
  
private:

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
  /** Compute shear velocities in case of relative slip with a given strength */
  inline void computeIndepShearVelocities(Real strength, UInt elem);
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

  /** Normal Strength */
  std::vector<Real> nor_strength;
  /** Shear Strength */
  std::vector<Real> shr_strength;
  /** Frictional strengh (when overlapping is prevented) */
  std::vector<Real> fric_strength;
  /** Cracking index just needed to better visualization of fracture process
      Standard value: 0 = outside the crack, 1 = in the cohesive zone,
      2 = inside the crack, 3 = inside the contact zone
      Other values can be defined arbitrarily at interface creation */
  std::vector<UInt> ind_crack;


  /** Critical normal opening */
  std::vector<Real> crit_nor_opening;
  /** Critical shear opening */
  std::vector<Real> crit_shr_opening;
  /** Maximum normal strength */
  std::vector<Real> max_nor_strength;
  /** Maximum shear strength */
  std::vector<Real> max_shr_strength;

  /** Overlapping tolerance (0 = no , 1 = yes) */
  bool allow_overlapping;
  // Associated contact law
  std::shared_ptr<ContactLaw> contact_law;

  /** Permanent access toward some fields registered in the DataRegister
      required to compute interface conditions */
  std::vector<Real> mu;
  std::vector<CrackProfile*> stresses;
  std::vector<CrackProfile*> velocities;
  std::vector<CrackProfile*> displacements;
  CrackProfile * intfc_trac;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                       

                    */
/* -------------------------------------------------------------------------- */

inline CohesiveLaw::~CohesiveLaw(){							  
										 
}							


#include "cohesive_law_inline_impl.cc"

#endif /* __COHESIVE_LAW__ */
