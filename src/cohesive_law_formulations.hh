/**
 * @file   cohesive_law__formulations.hh
 * @author Roxane Ferry <roxane.ferry@epfl.ch>
 * @date   Wed Feb 15 10:07:06 2023
 *
 * @brief  Struct describing the formulations of the cohesive laws
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
#include "cRacklet_common.hh"
#include "cohesive_law_all.hh"
/* -------------------------------------------------------------------------- */
struct CohesiveFormulation {
  CohesiveFormulation(){};
  virtual ~CohesiveFormulation(){};
  virtual inline Real getStrength(Real crit_nor_op, Real crit_shr_op, Real nor_op, Real shr_op, Real max_nor_str, 
  Real max_shr_str, Real res_nor_str, Real res_shr_str, Real &nor_str, Real &shr_str){
    Real aux;
    UInt id_crack;
    aux = sqrt((nor_op/crit_nor_op)*(nor_op/crit_nor_op)+
		(shr_op/crit_shr_op)*(shr_op/crit_shr_op));
    
    if ((aux>=1)||(nor_str==res_nor_str)||(shr_str==res_shr_str)) {

      bool in_contact = ((nor_str == res_nor_str)&&(cRacklet::is_negative((nor_op))));
      // the case of contact is handled by the associated ContactLaw
    
      if (!in_contact) { 
        id_crack = 2;
        nor_str = res_nor_str;
        shr_str= res_shr_str;
      }
    }

    else {
      id_crack = 1;
      nor_str = max_nor_str - (max_nor_str-res_nor_str)*aux;
      shr_str = max_shr_str - (max_shr_str-res_shr_str)*aux;
    }

    return id_crack;
  }
};

/* -------------------------------------------------------------------------- */
struct DualCohesiveFormulation : public CohesiveFormulation{
  DualCohesiveFormulation(){};
  virtual ~DualCohesiveFormulation(){};
  virtual inline Real getStrength(Real crit_nor_op, Real crit_shr_op, Real nor_op, Real shr_op, Real max_nor_str,
  Real max_shr_str, Real res_nor_str, Real res_shr_str, Real &nor_str, Real &shr_str){
    Real aux;
    UInt id_crack;
    Real nor_op_factor = CohesiveLawAll::nor_op_factor;
    Real shr_op_factor = CohesiveLawAll::shr_op_factor;
    Real nor_str_factor = CohesiveLawAll::nor_str_factor;
    Real shr_str_factor = CohesiveLawAll::shr_str_factor;

    aux = sqrt((nor_op/crit_nor_op)*(nor_op/crit_nor_op)+
		(shr_op/crit_shr_op)*(shr_op/crit_shr_op));

    if (aux/shr_op_factor <= 1) {
      id_crack = 1;
      nor_str = max_nor_str - (max_nor_str-res_nor_str*nor_str_factor)/(nor_op_factor)*aux;
      shr_str = max_shr_str - (max_shr_str-res_shr_str*shr_str_factor)/(shr_op_factor)*aux;
    }

    else if (aux < 1) {
      id_crack = 1;
      shr_str = (res_shr_str*(1-shr_str_factor)/(1-shr_op_factor)) * (aux - 1) + res_shr_str;
      nor_str = (res_nor_str*(1-nor_str_factor)/(1-nor_op_factor)) * (aux - 1) + res_nor_str;
    }

    else {
      bool in_contact = ((nor_str == res_nor_str)&&(cRacklet::is_negative((nor_op))));
      // the case of contact is handled by the associated ContactLaw
      if (!in_contact) { 
        id_crack = 2;
        nor_str = res_nor_str;
        shr_str = res_shr_str;
      }
    }

    return id_crack; 
  }
};
