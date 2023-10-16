/**
 * @file   cohesive_law_formulation_coulomb.hh
 * @author Roxane Ferry <roxane.ferry@epfl.ch>
 * @date   Tue Oct 3 11:16:47 2023
 *
 * @brief  Struct describing the cohesive coulomb formulations
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
#include "cohesive_law_coulomb.hh"
/* -------------------------------------------------------------------------- */
struct CohesiveCoulombFormulation {
  CohesiveCoulombFormulation(){};
  virtual ~CohesiveCoulombFormulation(){};
  virtual inline Real operator()(Real cf_s, Real cf_d, Real aux, Real shr_op, UInt i, UInt &id_crack){ 
    Real coeff;

    if (aux >= 1) {
      id_crack = 2;
      coeff = cf_d;
      std::cout << "Case D" << std::endl;
    }
    
    else {
    id_crack = 1;
    coeff = cf_s - (cf_s - cf_d) * aux;
    }
    return coeff;

  }
};

/* -------------------------------------------------------------------------- */
struct DualCohesiveCoulombFormulation : public CohesiveCoulombFormulation {
  DualCohesiveCoulombFormulation(){};
  virtual ~DualCohesiveCoulombFormulation(){};
  virtual inline Real operator()(Real cf_s, Real cf_d, Real aux, Real shr_op, UInt i, UInt &id_crack){
    std::vector<Real> cf_int = CohesiveLawCoulomb::cf_int;
    std::vector<Real> crit_int_op = CohesiveLawCoulomb::crit_int_opening;
    std::vector<Real> crit_shr_op = CohesiveLawCoulomb::crit_shr_opening;

    Real aux2 = shr_op/crit_int_op[i];
    Real coeff;

    if (shr_op <= crit_int_op[i]) {
      id_crack = 1;
      coeff = cf_s + ((cf_int[i] - cf_s) / (crit_int_op[i]))*shr_op; 
      std::cout << "Case B" << std::endl;
    }

    else if (aux < 1) {
      id_crack = 1;
      coeff = ((cf_int[i] - cf_d) / (crit_int_op[i] - crit_shr_op[i])) * (shr_op - crit_int_op[i]) + cf_int[i];
      std::cout << "Case C" << std::endl;
    }

    else {
      id_crack = 2;
      coeff = cf_d;
      std::cout << "Case D" << std::endl;
    }

          std::cout << "shr_op = " << shr_op << std::endl;
      std::cout << "crit_int_op[i] = " << crit_int_op[i] << std::endl;

    return coeff;
  }
};
