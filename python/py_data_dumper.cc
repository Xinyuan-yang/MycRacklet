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

    .def("initDumper",&DataDumper::initDumper,py::arg("filename"),py::arg("type"),py::arg("ratio_of_nele")=1.0,py::arg("stride")=1,py::arg("start")=0,py::arg("format")=_text,
	 "Create a Dumper to output a field of given type into filename, Eventually specify which amout (ratio) of the full data, the stride to appy while dumping and the first position to dump")

    .def("initVectorDumper",&DataDumper::initVectorDumper,py::arg("filename"), py::arg("type"), py::arg("dimension_to_dump"), py::arg("ratio_of_nele") = 1.0, py::arg("stride")=1, py::arg("start")=0, py::arg("format")=_text,
	 "Create a Dumper to output vectorial field (of size n_ele*dim). Options are the same than initDumper + one entry to specify which direction to dump")

    .def("initPointsDumper",[](DataDumper & self,const std::string filename, std::vector<DataFields> fields,
			       std::vector<UInt> points_to_dump, OutputFormat format=_text) {self.initPointsDumper(filename,fields,points_to_dump,format);},
	 "Create a PointsDumper to output a list of fields at precise interface position specified in points_to_dump. IMPORTANT: With binary mode, every data is written in the format of Real data type !")

    .def("initPointsDumper",[](DataDumper & self,const std::string filename, std::vector<UInt> points_to_dump, OutputFormat format=_text) {self.initPointsDumper(filename,points_to_dump,format);},
	 "Create a PointsDumper as above but using the standard output fields listed on top of this file")

    .def("initPointsDumper",[](DataDumper & self,const std::string filename, UInt start, UInt end, UInt nb_obs_points, OutputFormat format=_text) {self.initPointsDumper(filename,start,end,nb_obs_points,format);},
	 "Create a PointsDumper of standard fields along nb_obs_points between elements start and end")

    .def("initIntegratorsDumper",[](DataDumper & self,const std::string filename, std::vector<UInt> integration_domain,
			     std::vector<IntegratorTypes> inte_types=standard_energy_integrators,
			     std::vector<std::string> integrator_names=standard_energy_names,
				    OutputFormat format=_text) {self.initIntegratorsDumper(filename,integration_domain,inte_types,integrator_names,format);},
	 "Create a IntegratorsDumper to compute/output a list of integration types along a given integration domain. Each integration type gets an associated name used for registration in the DataRegister")
    
    .def("initIntegratorsDumper",[](DataDumper & self,const std::string filename, std::vector<Real> start_corner, std::vector<Real> end_corner, OutputFormat format=_text) {self.initIntegratorsDumper(filename,start_corner,end_corner,format);},
	 "Create a IntegratorsDumper with standard energies within a rectangular subset of the plane between start and en corner")
    
    .def("initIntegratorsDumper",[](DataDumper & self,const std::string filename, OutputFormat format=_text) {self.initIntegratorsDumper(filename,format);},
	 "Create a IntegratorsDumper with standard energies integrated over the entire plane")

    .def("initSurfingIntegratorsDumper",&DataDumper::initSurfingIntegratorsDumper,
	 py::arg("filename"),py::arg("integration_width"),py::arg("crack_start"),py::arg("crack_end"),py::arg("inte_types"),py::arg("integrator_names"),py::arg("format")=_text,
	 "Create a IntegratorsDumper which follows the propagating tip with an integration domain of a given width. Crack position is tracked between crack_start and crack_end")
    .def("dump",&DataDumper::dump,
	 py::arg("filename"),
	 "Generate an output of current model state to given file")
    .def("dumpAll",&DataDumper::dumpAll,
	 "Generate output from every created dumper at once");
}

} // namespace cracklet
