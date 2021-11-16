/**
 * @file   py_data_register.cc
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
#include <data_register.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cracklet {
/* -------------------------------------------------------------------------- */
  
void register_data_fields(py::module& mod) {
  py::enum_<DataFields>(mod,"DataFields")
    .value("_top_displacements",DataFields::_top_displacements)
    .value("_bottom_displacements",DataFields::_bottom_displacements)
    .value("_normal_displacement_jumps",DataFields::_normal_displacement_jumps)
    .value("_shear_displacement_jumps",DataFields::_shear_displacement_jumps)
    .value("_top_velocities",DataFields::_top_velocities)
    .value("_bottom_velocities",DataFields::_bottom_velocities)
    .value("_normal_velocity_jumps",DataFields::_normal_velocity_jumps)
    .value("_shear_velocity_jumps",DataFields::_shear_velocity_jumps)
    .value("_interface_tractions",DataFields::_interface_tractions)
    .value("_top_loading",DataFields::_top_loading)
    .value("_bottom_loading",DataFields::_bottom_loading)
    .value("_top_dynamic_stress",DataFields::_top_dynamic_stress)
    .value("_bottom_dynamic_stress",DataFields::_bottom_dynamic_stress)
    .value("_normal_strength",DataFields::_normal_strength)
    .value("_shear_strength",DataFields::_shear_strength)
    .value("_frictional_strength",DataFields::_frictional_strength)
    .value("_id_crack",DataFields::_id_crack)
    .value("_critical_normal_opening",DataFields::_critical_normal_opening)
    .value("_critical_shear_opening",DataFields::_critical_shear_opening)
    
    .value("_state_variable",DataFields::_state_variable)
    .value("_friction_coefficient",DataFields::_friction_coefficient)
    .value("_rands_D",DataFields::_rands_D)
    .value("_rands_f_0",DataFields::_rands_f_0)
    .value("_rands_a",DataFields::_rands_a)
    .value("_rands_b",DataFields::_rands_b)
    .value("_rands_v_star",DataFields::_rands_v_star)
    .value("_rands_phi_star",DataFields::_rands_phi_star)
    .export_values();

}

void register_integrator_types(py::module& mod) {
  py::enum_<IntegratorTypes>(mod,"IntegratorTypes")
    .value("_shear_fracture_energy",IntegratorTypes::_shear_fracture_energy)
    .value("_normal_fracture_energy",IntegratorTypes::_normal_fracture_energy)
    .value("_frictional_energy",IntegratorTypes::_frictional_energy)
    .value("_radiated_energy",IntegratorTypes::_radiated_energy)
    .export_values();
}  
    
void register_data_register(py::module& mod) {
  py::class_<DataRegister>(mod, "DataRegister")
    .def(py::init<>())
    .def_readwrite_static("restart_dir",&DataRegister::restart_dir)
    .def("readData",&DataRegister::readData)
    .def("registerParameterReal",&DataRegister::registerParameter<Real>)
    .def("registerParameterInt",&DataRegister::registerParameter<UInt>)
    .def("registerParameterString",&DataRegister::registerParameter<std::string>)
    .def("readInputFile",&DataRegister::readInputFile)
    .def("getParameterReal",&DataRegister::getParameter<Real>)
    .def("getParameterInt",&DataRegister::getParameter<UInt>)
    .def("getParameterString",&DataRegister::getParameter<std::string>)
    .def("registerComputer",&DataRegister::registerComputer)
    .def("getComputer",&DataRegister::getComputer)
    .def("getCrackTipPosition",&DataRegister::getCrackTipPosition,py::arg(),py::arg("z_pos")=0)
    .def("getTopVelocities",&DataRegister::getTopVelocities,py::return_value_policy::reference)
    .def("getBotVelocities",&DataRegister::getBotVelocities,py::return_value_policy::reference)
    .def("getShearVelocityJumps",&DataRegister::getShearVelocityJumps,py::return_value_policy::reference)
    .def("getNormalVelocityJumps",&DataRegister::getNormalVelocityJumps,py::return_value_policy::reference)
    .def("getTopDisplacements",&DataRegister::getTopDisplacements,py::return_value_policy::reference)
    .def("getBotDisplacements",&DataRegister::getBotDisplacements,py::return_value_policy::reference)
    .def("getShearDisplacementJumps",&DataRegister::getShearDisplacementJumps,py::return_value_policy::reference)
    .def("getNormalDisplacementJumps",&DataRegister::getNormalDisplacementJumps,py::return_value_policy::reference)
    .def("getInterfaceTractions",&DataRegister::getInterfaceTractions,py::return_value_policy::reference)
    .def("getTopDynamicStresses",&DataRegister::getTopDynamicStresses,py::return_value_policy::reference)
    .def("getBotDynamicStresses",&DataRegister::getBotDynamicStresses,py::return_value_policy::reference)
    .def("getTopLoading",&DataRegister::getTopLoading,py::return_value_policy::reference)
    .def("getBotLoading",&DataRegister::getBotLoading,py::return_value_policy::reference)
;
  
}

} // namespace cracklet
