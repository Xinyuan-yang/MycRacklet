/**
 * @file   py_data_dumper.cc
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
#include <data_dumper.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cracklet {
/* -------------------------------------------------------------------------- */

void register_output_format(py::module& mod) {
  py::enum_<OutputFormat>(mod,"OutputFormat")
    .value("_text",OutputFormat::_text)
    .value("_binary",OutputFormat::_binary);
  }
  
void register_data_dumper(py::module& mod) {
  
  py::class_<DataDumper>(mod, "DataDumper")
    .def(py::init<SpectralModel&>())

    .def("initDumper",&DataDumper::initDumper,py::arg("filename"),py::arg("type"),py::arg("ratio_of_nele")=1.0,py::arg("stride")=1,py::arg("start")=0,py::arg("format")=_text)

    .def("initVectorDumper",&DataDumper::initVectorDumper,py::arg("filename"), py::arg("type"), py::arg("dimension_to_dump"), py::arg("ratio_of_nele") = 1.0, py::arg("stride")=1, py::arg("start")=0, py::arg("format")=_text)

    .def("initPointsDumper",py::overload_cast<const std::string, std::vector<DataFields>, std::vector<UInt>, OutputFormat>(&DataDumper::initPointsDumper), py::arg("filename"), py::arg("fields"), py::arg("points_to_dump"), py::arg("format")=_text)

    .def("initPointsDumper",py::overload_cast<const std::string, std::vector<UInt>, OutputFormat>(&DataDumper::initPointsDumper),py::arg("filename"), py::arg("points_to_dump"), py::arg("format")=_text)

    .def("initPointsDumper",py::overload_cast<const std::string, UInt, UInt, UInt, OutputFormat>(&DataDumper::initPointsDumper),py::arg("filename"), py::arg("start"), py::arg("end"), py::arg("nb_obs_points"), py::arg("format")=_text)

    .def("initIntegratorsDumper",py::overload_cast<const std::string, std::vector<UInt>, std::vector<IntegratorTypes>, std::vector<std::string>,OutputFormat>(&DataDumper::initIntegratorsDumper),
	 py::arg("filename"), py::arg("integration_domain"), py::arg("inte_types")=standard_energy_integrators, py::arg("integrator_names")=standard_energy_names, py::arg("format")=_text)
    
    .def("initIntegratorsDumper",py::overload_cast<const std::string, std::vector<Real>, std::vector<Real>, OutputFormat>(&DataDumper::initIntegratorsDumper),
	 py::arg("filename"), py::arg("start_corner"), py::arg("end_corner"), py::arg("format")=_text)
    
    .def("initIntegratorsDumper",py::overload_cast<const std::string, OutputFormat>(&DataDumper::initIntegratorsDumper),
	 py::arg("filename"), py::arg("format")=_text)

    .def("initSurfingIntegratorsDumper",&DataDumper::initSurfingIntegratorsDumper)
    .def("dump",&DataDumper::dump)
    .def("dumpAll",&DataDumper::dumpAll);
}

} // namespace cracklet
