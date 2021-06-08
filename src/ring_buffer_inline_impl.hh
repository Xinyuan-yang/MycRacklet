/**
 * @file   ring_buffer_inline_impl.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Nov 18 09:42:13 2015
 *
 * @brief  Implementation of the inline functions of the RingBuffer class
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
template<typename T>
inline RingBuffer<T>::RingBuffer(){

  values = NULL;
  size = 0;
  n = 0;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline RingBuffer<T>::~RingBuffer(){

}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void RingBuffer<T>::init(T * field_val, UInt size, T default_val){

  this->size = size;
  this->values = field_val;
  n=0;
  T * it = values;
  for (UInt i = 0; i < size; ++i) {
    *it = default_val;
    ++it;
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void RingBuffer<T>::operator<<(T new_val){

  *(values+n%size) = new_val;
  ++n;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline T * RingBuffer<T>::current() const {
  return values+(n-1)%size;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline T * RingBuffer<T>::begin() const {
  return values;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline T * RingBuffer<T>::end() const {
  return values+std::min(size-1,n-1);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void RingBuffer<T>::printself(std::ostream & stream) const {

  T * it = values;
  for (UInt i = 0; i < size; ++i) {
    stream << " " << *it;
    ++it;
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline std::ostream & operator <<(std::ostream & stream, const RingBuffer<T> & _this)
{
  _this.printself(stream);
  return stream;
}
