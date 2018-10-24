/**
 * @file   interface_law.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Sun Jan  6 19:05:32 2013
 *
 * @brief  Abstract class representing the laws governing interface conditions
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
#ifndef __INTERFACE_LAW__
#define __INTERFACE_LAW__
/* -------------------------------------------------------------------------- */
#include "crack_profile.hh"
#include "data_register.hh"
#include <iostream>
/* -------------------------------------------------------------------------- */

class InterfaceLaw : public DataRegister {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  InterfaceLaw();
  virtual ~InterfaceLaw();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // initialize the conditions at the interface
  virtual void initInterfaceConditions() = 0;
  // solve the interface conditions in function of the opening profile and dynamic stress field
  virtual void updateInterfaceConditions() = 0;
  // Method used in restart framework in case of history dependant fracture law
  // pausing=true->generate restart files | pausing=false->restart simulation from existing files
  // If 3d simulation is restarted from 2d one, specify the number of elements along x (nele_2d=nele_x)
  virtual void restart(bool pausing=true, UInt nele_2d=0)=0;

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
inline InterfaceLaw::InterfaceLaw(){							
 										
}										
/* -------------------------------------------------------------------------- */
inline InterfaceLaw::~InterfaceLaw(){							  
										 
}										

#endif /* __INTERFACE_LAW__ */
 
  
