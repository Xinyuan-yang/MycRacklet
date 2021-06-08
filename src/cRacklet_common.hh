/**
 * @file   cRacklet_common.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Tue Feb  9 14:40:13 2016
 *
 * @brief  Common methods and types of cRacklet
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
#ifndef __CRACKLET_COMMON__
#define __CRACKLET_COMMON__
/* -------------------------------------------------------------------------- */
#include <exception>
#include <string>
#include <iostream>
#include <sstream>
#include <limits>
#include <cmath>
#include <memory>
/* -------------------------------------------------------------------------- */
//Generic Types
typedef unsigned int UInt;
typedef size_t Idx;
typedef double Real;
//Default stable time step parameter beta such as c_s*dt = beta*dx
static const double CS_DT_OVER_DX=0.2;
namespace cRacklet {
  
  // Error handling
  /* -------------------------------------------------------------------------- */
  static inline void error(const std::string & message) {
    std::string err = "!!! cRacklet Error !!! ";
    std::cerr << err+message << std::endl;
    throw;
  }
  
  /* -------------------------------------------------------------------------- */
  static inline void error(const std::stringstream & message) {
    error(message.str());
  }

  /* -------------------------------------------------------------------------- */
  static inline bool has_converged(const Real error) {
    return std::abs(error)<1e-13;
  }
  
  // Compare is two values are identical up to machine precision
  /* -------------------------------------------------------------------------- */
  template<typename T>
  static inline bool is_equal(const T value1, const T value2) {
    
    return std::abs((value1-value2)/value1) < std::numeric_limits<T>::epsilon();
  }
  
  // Check if a value is negative up to machine precision
  /* -------------------------------------------------------------------------- */
  template<typename T>
  static inline bool is_negative(const T value) {

    return ((std::abs(value) > std::numeric_limits<T>::epsilon())&&(value<0.0)); 
  }
}

#endif /* __CRACKLET_COMMON__ */
