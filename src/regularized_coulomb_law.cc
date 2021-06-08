/**
 * @file   regularized_coulomb_law.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Mar 18 21:56:18 2013
 *
 * @brief  Implementation of the RegularizedCoulombLaw class
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
#include "regularized_coulomb_law.hh"
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
/* -------------------------------------------------------------------------- */
void RegularizedCoulombLaw::computeFricStrength(Real & norm_comp_stress, Real & strength, UInt i, UInt it) {


  if (it == 1) {

    sigma_np1[i] = norm_comp_stress;

  }

  else if ((contact_history[i]+1) == it) {

   sigma_np1[i] = (sigma_np1[i] + dt_tstar*norm_comp_stress) / (1 + dt_tstar);

   //norm_comp_stress = sigma_np1[i];

  }

  else {

    sigma_np1[i] = dt_tstar*norm_comp_stress/(1 + dt_tstar);

    //norm_comp_stress = sigma_np1[i];
  }

  strength = cf*fabs(sigma_np1[i]);
  contact_history[i] = it;

}

/* -------------------------------------------------------------------------- */
void RegularizedCoulombLaw::restart(bool pausing, UInt nele_2d) {

  restartData(sigma_np1, "restart_sigma_np1.cra",pausing, nele_2d);
  restartData(contact_history, "restart_contact_history.cra",pausing, nele_2d);
}

