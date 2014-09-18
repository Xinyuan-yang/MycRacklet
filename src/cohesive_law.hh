/**
 * @file   cohesive_law.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Sun Jan  6 19:35:47 2013
 *
 * @brief  Class describing a cohesive law
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
#ifndef __COHESIVE_LAW__
#define __COHESIVE_LAW__
/* -------------------------------------------------------------------------- */
#include "fracture_law.hh"
#include <vector>
#include <math.h>
#define PI 3.141592653589 
/* -------------------------------------------------------------------------- */

class CohesiveLaw : public FractureLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  CohesiveLaw(double crit_n_open, double crit_shr_open, double max_nor_strght, 
	      double max_shr_strght);

  virtual ~CohesiveLaw();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // update the strength of the material in function of the opening profile 
  void updateFractureLaw(std::vector<double> & nor_strength, std::vector<double> & shr_strength,
			 std::vector<unsigned int> & ind_crack, CrackProfile & nor_opening, 
			 CrackProfile & shr_opening);
  // dump current cohesive law in a given ofstream
  void printSelf(std::ofstream & parameters_file, std::ofstream & summary);

private:
 
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
   // Critical normal opening
  double crit_nor_opening;
  // Critical shear opening
  double crit_shr_opening;
  // Maximum normal strength
  double max_nor_strength;
  // Maximum shear strength
  double max_shr_strength;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline CohesiveLaw::CohesiveLaw(double crit_n_open, double crit_shr_open, double max_nor_strght,
				double max_shr_strght) : FractureLaw() {

  crit_nor_opening = crit_n_open;
  crit_shr_opening = crit_shr_open;
  max_nor_strength = max_nor_strght;
  max_shr_strength = max_shr_strght;
										
}										
/* -------------------------------------------------------------------------- */
inline CohesiveLaw::~CohesiveLaw(){							  
										 
}					


#endif /* __COHESIVE_LAW__ */
