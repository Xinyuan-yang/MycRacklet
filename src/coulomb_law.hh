/**
 * @file   coulomb_law.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Jan  7 13:51:56 2013
 *
 * @brief  Class describing coulomb contact law
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
#ifndef __COULOMB_LAW__
#define __COULOMB_LAW__
/* -------------------------------------------------------------------------- */
#include "contact_law.hh"
/* -------------------------------------------------------------------------- */

/**
 * @class  CoulombLaw coulomb_law.hh
 *
 * Class describing coulomb contact law
 *
*/
class CoulombLaw : public ContactLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  /** Constructor
      @param coef: (Real) coefficient of friction
   */
  CoulombLaw(Real coef);
  /// Default Destructor
  virtual ~CoulombLaw();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /** Compute the frictionnal strength with the normal compressive stress
      @param norm_comp_stress: (Real) value of the compressive stress
      @param strength: (Real) strength
      @param  i and it are unused
   */
  void computeFricStrength(Real & norm_comp_stress, Real & strength, UInt i, UInt it);
  /// Method used in restart framework but no history-dependant variable within this law
  void restart(bool pausing=false, UInt nele_2d=0){};
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
 
  // Coulomb coefficient of friction
  Real cf;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline CoulombLaw::CoulombLaw(Real coef) : ContactLaw(){							

  cf = coef;

  out_summary << "/* -------------------------------------------------------------------------- */ "; 
  out_summary << std::endl;
  out_summary << " CONTACT LAW VARIABLES " << std::endl;
  out_summary << "* Type of contact law: Coulomb law" << std::endl;		        
  out_summary << "* Coefficient of friction: " << cf << std::endl;
  out_summary << std::endl;	

  out_parameters << cf << " "; 
}										
/* -------------------------------------------------------------------------- */
inline CoulombLaw::~CoulombLaw(){							  
										 
}	


#endif /* __COULOMB_LAW__ */
