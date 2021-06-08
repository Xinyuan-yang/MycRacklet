/**
 * @file   regularized_coulomb_law.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Mar 18 21:56:18 2013
 *
 * @brief  Class computing regularized coulomb friction law
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
#ifndef __REGULARIZED_COULOMB_LAW__
#define __REGULARIZED_COULOMB_LAW__
/* -------------------------------------------------------------------------- */
#include "contact_law.hh"
#include <vector>
/* -------------------------------------------------------------------------- */

/**
 * @class  RegularizedCoulombLaw regularized_coulomb_law.hh
 *
 * Class computing regularized coulomb friction law
 *
*/
class RegularizedCoulombLaw : public ContactLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  /** Constructor
      @param coef: (Real) coefficient of friction
      @param ratio: (Real) regularization ratio (dt/t*)
      @param n_ele: (UInt) number of interface points
   */
  RegularizedCoulombLaw(Real coef, Real ratio, UInt n_ele);
  /// Default Destructor
  virtual ~RegularizedCoulombLaw();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /** Compute the frictionnal strength with the normal compressive stress
      @param norm_comp_stress: (Real) value of the compressive stress
      @param strength: (Real) strength
      @param  i: (UInt) element
      @param it: (UInt) iteration number
   */
  void computeFricStrength(Real & norm_comp_stress, Real & strength, UInt i, UInt it); 
  /// Method used in restart framework of contact history
  /// pausing=true->generate restart files | pausing=false->restart simulation from existing files
  /// If 3d simulation is restarted from 2d one, specify the number of elements along x (nele_2d=nele_x)
  void restart(bool pausing=false, UInt nele_2d=0);
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
  // Regularization ratio (dt/t*)
  Real dt_tstar;
  // Regularized contact pressure history
  std::vector<Real> sigma_np1;
  // Save the # of time step when element i are under contact
  std::vector<UInt> contact_history;
  
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline RegularizedCoulombLaw::RegularizedCoulombLaw(Real coef, Real ratio, UInt n_ele) : ContactLaw(){

  cf = coef;
  dt_tstar = ratio;
  sigma_np1.resize(n_ele);
  contact_history.resize(n_ele);

  out_summary << "/* -------------------------------------------------------------------------- */ "; 
  out_summary << std::endl;
  out_summary << " CONTACT LAW VARIABLES " << std::endl;
  out_summary << "* Type of contact law: Coulomb law with normal pressure regularization" << std::endl;	        
  out_summary << "* Coefficient of friction: " << cf << std::endl;
  out_summary << "* Regularization time scale dt/t*:  " << dt_tstar << std::endl; 
  out_summary << std::endl;	

  out_parameters << cf << " ";

  
}										
/* -------------------------------------------------------------------------- */
inline RegularizedCoulombLaw::~RegularizedCoulombLaw(){
										 
}	


#endif /* __REGULARIZED_COULOMB_LAW__ */
