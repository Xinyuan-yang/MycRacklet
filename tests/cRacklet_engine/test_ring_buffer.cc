/**
 * @file   test_ring_buffer.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Nov 19 09:30:44 2015
 *
 * @brief  Testing ring_buffer objects
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
#include "ring_buffer.hh"
#include <iostream>
/* -------------------------------------------------------------------------- */

int main() {

  UInt buffer_size = 20;
  UInt nb_steps = 50;

  RingBuffer<int> buffer;

  int * val = new int[buffer_size];

  buffer.init(val, buffer_size,-1);

  int * it;
  int * current;
  int * begin;
  int * end;

  for (UInt t = 0; t < nb_steps; ++t) {

    buffer << t;
    
    std::cout << buffer << std::endl;

    current = buffer.current()+1;
    begin = buffer.begin();
    end = buffer.end()+1;

    for(it = current; it!=end; ++it) {

      std::cout << " "<< *it;
    }

    std::cout << " | ";
    
    for(it = begin; it!=current; ++it){

      std::cout << " "<< *it;
    }
    
    std::cout << std::endl << std::endl;

  }
  
  delete [] val;

  return 0;
}
