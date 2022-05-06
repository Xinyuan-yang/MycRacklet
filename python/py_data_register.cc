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
    .value("_maximum_normal_strength",DataFields::_maximum_normal_strength)
    .value("_maximum_shear_strength",DataFields::_maximum_shear_strength)
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

void register_data_types(py::module& mod) {
  py::class_<DataTypes>(mod, "DataTypes", py::buffer_protocol())
    .def(py::init<>())
    .def("testMe",[](DataTypes&datatype){
	UInt res;
	if(datatype.vec_double != NULL)
	  res=1;
	else if(datatype.vec_uint != NULL)
	  res=2;
	else if(datatype.crack_prof != NULL)
	  res=3;
	else
	  res=0;
	return res;
      })
    .def_buffer([](DataTypes&datatype) -> py::buffer_info {
	if(datatype.vec_double != NULL)
	  return py::buffer_info(datatype.vec_double->data(),
				 sizeof(Real),
				 py::format_descriptor<Real>::format(),
				 1,
				 {datatype.vec_double->size()},
				 {sizeof(Real)}
				 );
	else if(datatype.vec_uint != NULL)
	  return py::buffer_info(datatype.vec_uint->data(),
				 sizeof(UInt),
				 py::format_descriptor<UInt>::format(),
				 1,
				 {datatype.vec_uint->size()},
				 {sizeof(UInt)}
				 );
	else if(datatype.crack_prof != NULL) {
	  std::vector<UInt> nele = datatype.crack_prof->shape();
	  UInt dim = datatype.crack_prof->size()/(nele[0]*nele[1]);
	  return py::buffer_info(datatype.crack_prof->getValues().data(),
				 sizeof(Real),
				 py::format_descriptor<Real>::format(),
				 3,
				 {nele[0],nele[1],dim},
				 {sizeof(Real)*dim,sizeof(Real)*dim*nele[0],sizeof(Real)}
				 );
	}
	else
	  return py::buffer_info();
    });
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
    .def("readData",&DataRegister::readData,
	 "Access to a given field returned in the type DataTypes")
    .def("registerParameterReal",&DataRegister::registerParameter<Real>,
	 "Register a global simulation parameter, whose type is a float")
    .def("registerParameterInt",&DataRegister::registerParameter<UInt>,
	 "Register a global simulation parameter, whose type is an integer")
    .def("registerParameterString",&DataRegister::registerParameter<std::string>,
	 "Register a global simulation parameter, whose type is a string")
    .def("readInputFile",&DataRegister::readInputFile,
	 "Register a set of global parameters from input file")
    .def("getParameterReal",&DataRegister::getParameter<Real>,
	 "Get the value of a simulation parameter, whose type is a float")
    .def("getParameterInt",&DataRegister::getParameter<UInt>,
	 "Get the value of a simulation parameter, whose type is an integer")
    .def("getParameterString",&DataRegister::getParameter<std::string>,
	 "Get the value of a simulation parameter, whose type is a string")
    .def("registerComputer",&DataRegister::registerComputer,
	 "Register a computer object with a given name")
    .def("getComputer",&DataRegister::getComputer,
	 "Access a registered computer object")
    .def("getCrackTipPosition",&DataRegister::getCrackTipPosition,py::arg(),py::arg("z_pos")=0,
	 "Method returning current crack position searched between x_start and x_end (looking for the first point with state == 2)")
    .def("getCohesiveTipPosition",&DataRegister::getCohesiveTipPosition,py::arg(),py::arg("z_pos")=0,
	 "Method returning current cohesive tip position searched between x_start and x_end (looking for the last point with state == 1)")
    .def("getTopVelocities",&DataRegister::getTopVelocities,py::return_value_policy::reference,
	 "Direct access to top velocities")
    .def("getBotVelocities",&DataRegister::getBotVelocities,py::return_value_policy::reference,
	 "Direct access to bottom velocities")
    .def("getShearVelocityJumps",&DataRegister::getShearVelocityJumps,py::return_value_policy::reference,
	 "Direct access to shear velocity jumps")
    .def("getNormalVelocityJumps",&DataRegister::getNormalVelocityJumps,py::return_value_policy::reference,
	 "Direct access to normal velocity jumps")
    .def("getTopDisplacements",&DataRegister::getTopDisplacements,py::return_value_policy::reference,
	 "Direct access to top displacements")
    .def("getBotDisplacements",&DataRegister::getBotDisplacements,py::return_value_policy::reference,
	 "Direct access to bottom displacements")
    .def("getShearDisplacementJumps",&DataRegister::getShearDisplacementJumps,py::return_value_policy::reference,
	 "Direct access to shear displacements jumps")
    .def("getNormalDisplacementJumps",&DataRegister::getNormalDisplacementJumps,py::return_value_policy::reference,
	 "Direct access to normal displacements jumps")
    .def("getInterfaceTractions",&DataRegister::getInterfaceTractions,py::return_value_policy::reference,
	 "Direct access to interface tractions")
    .def("getTopDynamicStresses",&DataRegister::getTopDynamicStresses,py::return_value_policy::reference,
	 "Direct access to top dynamic stresses")
    .def("getBotDynamicStresses",&DataRegister::getBotDynamicStresses,py::return_value_policy::reference,
	 "Direct access to bottom dynamic stresses")
    .def("getTopLoading",&DataRegister::getTopLoading,py::return_value_policy::reference,
	 "Direct access to top loading")
    .def("getBotLoading",&DataRegister::getBotLoading,py::return_value_policy::reference,
	 "Direct access to bottom loading")
;
  
}

} // namespace cracklet
