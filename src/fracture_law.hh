/**
 * @file   material.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Sun Jan  6 19:05:32 2013
 *
 * @brief  Mother class dealing with fracure's law of crack interface
 *
 * @section LICENSE
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
#ifndef __FRACTURE_LAW__
#define __FRACTURE_LAW__
/* -------------------------------------------------------------------------- */
#include "crack_profile.hh"
#include <iostream>
/* -------------------------------------------------------------------------- */

class FractureLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  FractureLaw();
  virtual ~FractureLaw();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // update the strength of the material in function of the opening profile
  virtual void updateFractureLaw(std::vector<Real> & nor_strength, std::vector<Real> & shr_strength,
				 std::vector<UInt> & ind_crack, CrackProfile & nor_opening, 
				 CrackProfile & shr_opening) = 0;
  // Method used in restart framework in case of history dependant fracture law
  // pausing=true->generate restart files | pausing=false->restart simulation from existing files
  // If 3d simulation is restarted from 2d one, specify the number of elements along x (nele_2d=nele_x)
  virtual void restart(bool pausing=true, UInt nele_2d=0)=0;
  // dump fracture law paramters in a given ofstream
  virtual void printSelf(std::ofstream & parameters_file, std::ofstream & summary) = 0;
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
 
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline FractureLaw::FractureLaw(){							
 										
}										
/* -------------------------------------------------------------------------- */
inline FractureLaw::~FractureLaw(){							  
										 
}										

#endif /* __FRACTURE_LAW__ */
 
  
