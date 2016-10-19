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

/* -------------------------------------------------------------------------- */
#ifndef __COULOMB_LAW__
#define __COULOMB_LAW__
/* -------------------------------------------------------------------------- */
#include "contact_law.hh"
/* -------------------------------------------------------------------------- */

class CoulombLaw : public ContactLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  CoulombLaw(Real coef);
  virtual ~CoulombLaw();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // Compute the frictionnal strength with the normal compressive stress
  void computeFricStrength(Real & norm_comp_stress, Real & strength, UInt i, UInt it);
  // dump current contact law in a given ofstream
  void printSelf(std::ofstream & parameters_file, std::ofstream & summary);
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
										
}										
/* -------------------------------------------------------------------------- */
inline CoulombLaw::~CoulombLaw(){							  
										 
}	


#endif /* __COULOMB_LAW__ */
