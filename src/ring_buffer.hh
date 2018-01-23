/**
 * @file   ring_buffer.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Nov 18 09:42:13 2015
 *
 * @brief  Object handling First-In First-Out (FIFO) storage needed for convolution
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
 */
/* -------------------------------------------------------------------------- */
#ifndef __RING_BUFFER__
#define __RING_BUFFER__
/* -------------------------------------------------------------------------- */
#include <vector>
#include <cstddef>
#include <iostream>
#include "cRacklet_common.hh"
/* -------------------------------------------------------------------------- */

template<typename T>
class RingBuffer {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  RingBuffer();
  virtual ~RingBuffer();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // Initialize the buffer on an array of type T (already allocated) and length size
  // filled with a default value
  inline void init(T * field_val, UInt size, T default_val);
  // Insert a new value in the buffer
  inline void operator<<(T new_val);
  // Get the number of inserted element since buffer initialization
  inline UInt getStep(){return n;}
  // Reset the number of insert element since buffer initialization.
  // Mostly used when restarting a buffer
  inline void resetStep(UInt step){this->n=step;}
  // Address of the last value inserted in the buffer
  inline T * current() const;
  // Address of the first element in memory
  inline T * begin() const;
  // Address of the last element in memory
  inline T * end() const;
  // Output values contained in the buffer with same configuration than in memory
  inline void printself(std::ostream & stream) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  
  // Buffer size
  UInt size;
  // Pointer to the array of T allocated for FIFO operations
  T * values;
  // Number of inserted element since buffer initialization
  UInt n;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "ring_buffer_inline_impl.cc"

#endif /* __RING_BUFFER__ */
