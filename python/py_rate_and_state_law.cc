/**
 * @file   py_rate_and_state_law.cc
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
#include <rate_and_state_law.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cracklet {
/* -------------------------------------------------------------------------- */
    
void register_rate_and_state_law(py::module& mod) {
  py::class_<RateAndStateLaw,InterfaceLaw>(mod, "RateAndStateLaw")
    .def(py::init<>())
    // Init the friction laws
    .def("initStandardFormulation",&RateAndStateLaw::initStandardFormulation)
    .def("initVelocityWeakeningFormulation",&RateAndStateLaw::initVelocityWeakeningFormulation)
    .def("initRegularizedFormulation",&RateAndStateLaw::initRegularizedFormulation)
    // Init the state laws
    .def("initStateEvolution",&RateAndStateLaw::initStateEvolution)
    .def("initRegularizedStateEvolution",&RateAndStateLaw::initRegularizedStateEvolution)
    .def("initStateEvolution",&RateAndStateLaw::initStateEvolution)

    .def("setVelocityPredictor",&RateAndStateLaw::setVelocityPredictor)
    .def("initInterfaceConditions",&RateAndStateLaw::initInterfaceConditions)
    .def("updateInterfaceConditions",&RateAndStateLaw::updateInterfaceConditions)
    .def("perturbState",&RateAndStateLaw::perturbState)
    .def("insertPerturbationPatch",&RateAndStateLaw::insertPerturbationPatch)
    .def("insertGaussianPerturbation",&RateAndStateLaw::insertGaussianPerturbation)
    .def("restart",&RateAndStateLaw::restart);
}
} // namespace cracklet
