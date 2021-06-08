/**
 * @file   py_cohesive_law_viscoelastic.cc
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
#include <cohesive_law_viscoelastic.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cracklet {
/* -------------------------------------------------------------------------- */
    
void register_cohesive_law_viscoelastic(py::module& mod) {
  py::class_<CohesiveLawViscoelastic,CohesiveLaw>(mod, "CohesiveLawViscoelastic")
    .def(py::init<>())
    .def("initLinearFormulation",&CohesiveLawViscoelastic::initLinearFormulation)
    .def("initQuadraticFormulation",&CohesiveLawViscoelastic::initQuadraticFormulation)
    .def("initPowerLawFormulation",&CohesiveLawViscoelastic::initPowerLawFormulation)
    .def("initInterfaceConditions",&CohesiveLawViscoelastic::initInterfaceConditions)
    .def("updateInterfaceConditions",&CohesiveLawViscoelastic::updateInterfaceConditions)
    .def("restart",&CohesiveLawViscoelastic::restart);      
}

} // namespace cracklet
