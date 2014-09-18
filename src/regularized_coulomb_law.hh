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
#ifndef __REGULARIZED_COULOMB_LAW__
#define __REGULARIZED_COULOMB_LAW__
/* -------------------------------------------------------------------------- */
#include "contact_law.hh"
#include <vector>
/* -------------------------------------------------------------------------- */

class RegularizedCoulombLaw : public ContactLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  RegularizedCoulombLaw(double coef, double ratio, int n_ele);
  virtual ~RegularizedCoulombLaw();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // Compute the frictionnal strength with the normal compressive stress
  void computeFricStrength(double & norm_comp_stress, double & strength, int i, int it); 
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
  double cf;
  // Regularization ratio (dt/t*)
  double dt_tstar;
  // Regularized contact pressure history
  std::vector<double> sigma_np1;
  // Save the # of time step when element i are under contact
  std::vector<int> contact_history;
  
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline RegularizedCoulombLaw::RegularizedCoulombLaw(double coef, double ratio, int n_ele) : ContactLaw(){							

  cf = coef;
  dt_tstar = ratio;
  sigma_np1.resize(n_ele);
  contact_history.resize(n_ele);
										
}										
/* -------------------------------------------------------------------------- */
inline RegularizedCoulombLaw::~RegularizedCoulombLaw(){							  
										 
}	


#endif /* __REGULARIZED_COULOMB_LAW__ */
