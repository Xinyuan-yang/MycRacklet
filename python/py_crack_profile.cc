/**
 * @file   py_crack_profile.cc
 * @author Thibault Roch <thibault.roch@epfl.ch>
 *
 * @section LICENSE
 *
 * cRacklet - A spectral boundary integral method for interface fracture simulation
 * Copyright (©) 2012 - 2013 Fabian Barras
 *               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 *               Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include <crack_profile.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cracklet {
/* -------------------------------------------------------------------------- */
    
  void register_crack_profile(py::module& mod) {
    py::class_<CrackProfile>(mod, "CrackProfile")
      .def(py::init<>())
      .def(py::init<std::vector<UInt>,UInt>())
      .def("__call__",[](CrackProfile&crack, int i){return crack[i];})
      .def("__call__",[](CrackProfile&crack, py::array_t<UInt> idx){
			py::buffer_info buf = idx.request();
			auto profile = py::array_t<Real>(buf.size);
			py::buffer_info pro = profile.request();
			
			Real *ptr1 = static_cast<Real *>(pro.ptr);
			UInt *ptr2 = static_cast<UInt *>(buf.ptr);
			
			for (size_t i = 0; i < buf.shape[0]; i++){
			  ptr1[i] = crack[ptr2[i]];
			}
			return profile;});      
  }
  
} // namespace cracklet
