#include <data_register.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cRacklet {
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
    .def("getCrackTipPosition",&DataRegister::getCrackTipPosition)
    .def("getTopVelocities",&DataRegister::getTopVelocities)
    .def("getBotVelocities",&DataRegister::getBotVelocities)
    .def("getShearVelocityJumps",&DataRegister::getShearVelocityJumps)
    .def("getNormalVelocityJumps",&DataRegister::getNormalVelocityJumps)
    .def("getTopDisplacements",&DataRegister::getTopDisplacements)
    .def("getBotDisplacements",&DataRegister::getBotDisplacements)
    .def("getShearDisplacementJumps",&DataRegister::getShearDisplacementJumps)
    .def("getNormalDisplacementJumps",&DataRegister::getNormalDisplacementJumps)
    .def("getInterfaceTractions",&DataRegister::getInterfaceTractions);
  
}

} // namespace cRacklet
