/**
 * @file   cohesive_law_all.hh
 * @author Roxane Ferry <roxane.ferry@epfl.ch>
 * @date   Wed Feb 15 14:36:15 2023
 *
 * @brief  Class describing a coupled slip/opening weakening interface
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
#ifndef __COHESIVE_LAW_ALL__
#define __COHESIVE_LAW_ALL__
/* -------------------------------------------------------------------------- */
#include "cohesive_law.hh"
/* -------------------------------------------------------------------------- */

class CohesiveLawAll : public CohesiveLaw {

  #include "cohesive_law_formulations.hh"

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:    
  CohesiveLawAll()
    :CohesiveLaw() {}

  virtual ~CohesiveLawAll() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /** regular slip weakening law*/
  void initRegularFormulation();
  /** dual-scale slip weakening law*/
  void initDualFormulation();

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
  //inline void computeIndepNormalVelocities(UInt ix, UInt iz);
  /** Compute shear velocities with a given shear strength */
  //inline void computeShearVelocities(Real strength, UInt elem);
  /** Compute velocities in the case of contact at crack type */
  //inline void computeContactVelocities(UInt ix, UInt iz);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // Abstract object representing the associated viscoelastic formulation
  std::shared_ptr<CohesiveFormulation> formulation;
};

#endif /* __COHESIVE_LAW_ALL__ */
