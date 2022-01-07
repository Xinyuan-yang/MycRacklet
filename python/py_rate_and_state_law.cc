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
    .def("initStandardFormulation",&RateAndStateLaw::initStandardFormulation,
	 "Init a given R&S formulation and state evolution standard rate and state formulation")
    .def("initVelocityWeakeningFormulation",&RateAndStateLaw::initVelocityWeakeningFormulation,
	 "Init formulation with a unique velocity-weakening branch")
    .def("initRegularizedFormulation",&RateAndStateLaw::initRegularizedFormulation,
	 py::arg("v0"),
	 "Init rate and state formulation with a regularized stick-to-slip transition")
    // Init the state laws
    .def("initStateEvolution",&RateAndStateLaw::initStateEvolution,"Init the state evolution law")
    .def("initRegularizedStateEvolution",&RateAndStateLaw::initRegularizedStateEvolution,"Init Regularized state evolution")
    .def("initStateEvolution",&RateAndStateLaw::initStateEvolution,"Init state evolution")

    .def("setVelocityPredictor",&RateAndStateLaw::setVelocityPredictor,
	 py::arg("v_0_pred"),
	 "Define the velocity prediction for each component (used before searching the initial steady state)")
    .def("initInterfaceConditions",&RateAndStateLaw::initInterfaceConditions,
	 "Initialize interface fields")
    .def("updateInterfaceConditions",&RateAndStateLaw::updateInterfaceConditions,
	 "Update the strength of the material in function of the opening profile")
    // Compute the next velocity for imposed vBC
    .def("computeNextAverageVelocity",&RateAndStateLaw::computeNextAverageVelocity,
	 "Compute the interface conditions but do not update them! (Used in pseudo-velocity driven systems) and return the average velocity of the next (fictious) step.")
    .def("perturbState",[](RateAndStateLaw & self, Real epsilon, Real k) {self.perturbState(epsilon,k);},
	 "perturb the state variable by adding a sinusoidale perturbation")
    .def("perturbState",[](RateAndStateLaw & self, std::vector<Real> perturbation) {self.perturbState(perturbation);},
	 "perturbe the state variable by adding a vector to it")
    .def("insertPerturbationPatch",&RateAndStateLaw::insertPerturbationPatch,
	 py::arg("patch_limits"),py::arg("new_rate"),
	 "Insert gaussian perturbation patch in the velocity field")
    .def("insertGaussianPerturbation",&RateAndStateLaw::insertGaussianPerturbation,
	 py::arg("std_dev"),py::arg("amplitude"),
	 "Insert gaussian perturbation in the velocity field")
    .def("insertSkewedPerturbation",[](RateAndStateLaw & self, Real std_dev, Real amplitude, Real alpha, Real rel_loc = 0.5) {self.insertSkewedPerturbation(std_dev,amplitude,alpha,rel_loc); },
	 "Insert a skewed perturbation in the velocity field. The maximum of the perturbation is given by the amplitude")
    .def("insertSkewedPerturbation",[](RateAndStateLaw & self, std::vector<Real> std_dev, std::vector<Real> amplitude, std::vector<Real> alpha, std::vector<Real> rel_loc) {self.insertSkewedPerturbation(std_dev,amplitude,alpha,rel_loc);},
	 "Insert a skewed perturbation in the velocity field. The maximum of the perturbation is given by the amplitude")
    .def("restart",&RateAndStateLaw::restart,
	 "Method used in restart framework");
}
} // namespace cracklet
