/**
 * @file   py_simulation_driver.cc
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
#include <simulation_driver.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cracklet {
  /* -------------------------------------------------------------------------- */
  
  void register_load_control_type(py::module& mod) {
    py::enum_<LoadControlType>(mod,"LoadControlType")
      .value("_space_control",LoadControlType::_space_control)
      .value("_time_control",LoadControlType::_time_control)
      .export_values();
  }
  
  void register_simulation_driver(py::module& mod) {
    py::class_<SimulationDriver,DataRegister>(mod, "SimulationDriver")
      .def(py::init<SpectralModel&, Real>(),py::arg("model_to_drive"),py::arg("beta")=0.0,"Constructor for the Simulation Driver")
      .def(py::init<SpectralModel&, Real, Real>())
      .def("initConstantLoading",&SimulationDriver::initConstantLoading,py::arg("cst_loading"),py::arg("psi"),py::arg("phi"),"Init a simulation with constant loading conditions")
      .def("initLoadingFromFile",&SimulationDriver::initLoadingFromFile,py::arg("loading_file"),py::arg("load_control"),py::arg("initial_loading"),py::arg("psi"),py::arg("phi"),"Init a simulation with evolving loading conditions following loading_file")
      .def("initConstantSpeed",&SimulationDriver::initConstantSpeed,py::arg("initial_loading"),py::arg("psi"),py::arg("phi"),py::arg("average_max_stress"),py::arg("spont_crack_length")=0.0,py::arg("load_control")=LoadControlType::_time_control,py::arg("load_upper_bound")=0.9,py::arg("griffith_length")=0.,"Initialization before tailoring loading conditions to fix crack speed")
      .def("solveStep",&SimulationDriver::solveStep,"Solve one time step of the simulation")
      .def("writeLoading",&SimulationDriver::writeLoading,py::arg("load_file"),"Print tailored loading conditions to file")
      .def("launchCrack",&SimulationDriver::launchCrack,py::arg("crack_start"),py::arg("launched_size"),py::arg("v_init"),py::arg("one_side_propagation")=true,"Launching artificially a through crack from position crack_start up to a given launched_size by artificially growing it at speed v_init*cs");
  }
  
} // namespace cracklet
