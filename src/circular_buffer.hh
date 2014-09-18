/**
 * @file   circular_buffer.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Fri Nov 10 10:10:41 2012
 *
 * @brief   "Librarian" class storing variables history and computing convolution
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
#ifndef __CIRC_BUFFER__
#define __CIRC_BUFFER__
/* -------------------------------------------------------------------------- */
#include <vector>
/* -------------------------------------------------------------------------- */

template<class Bed>
class CircularBuffer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  CircularBuffer();
  CircularBuffer(std::vector<int> cut);
  virtual ~CircularBuffer();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  
  // Compute the convolution at t=it between kernel K and this buffer
  void convolution(std::vector<std::vector<double> > & K, std::vector<std::complex<double> > & res, int it);
  // store the i-th time step Bed in the buffer
  void store(Bed * datas, int it);
  // resize the Buffer
  void resize(std::vector<int> cut);
  // resize the Buffer and his members of type Bed
  void resize(int buf_s, int bed_s);
  // access to the it-th time step storage of the jth element
  void recall(std::vector<std::complex<double> > & res, int it);

private:
 
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  
private:
  
  // Buffer core of type Bed
  std::vector<std::vector<Bed> > box;
  // number of elements ( = nb of covolution to compute)
  int nele;


};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "circular_buffer_impl.cc"

#endif /* __CIRC_BUFFER__ */
