/**
 * @file   py_spectral_model.cc
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
#include <spectral_model.hh>
#include <interface_law.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cracklet {
/* -------------------------------------------------------------------------- */

void register_directions(py::module& mod) {
  py::enum_<SpatialDirection>(mod,"SpatialDirection")
    .value("_x",SpatialDirection::_x)
    .value("_y",SpatialDirection::_y)
    .value("_z",SpatialDirection::_z);
  }
  
void register_spectral_model(py::module& mod) {
  py::class_<SpectralModel,DataRegister>(mod, "SpectralModel")
    .def(py::init<>())
    .def(py::init<UInt, UInt, Real,
	 Real, Real, Real, UInt,
	 const std::string, const std::string>(),
	 py::arg("nele"), py::arg("nb_time_steps"),
	 py::arg("dom_size"), py::arg("nu"),
	 py::arg("E"), py::arg("cs"), py::arg("tcut"),
	 py::arg("simulation_summary"), py::arg("output_dir")="./",
	 "Constructor for 1D (line) interface between similar bodies")
    .def(py::init<UInt, UInt, Real,
	 Real, Real, Real, Real, Real, Real, UInt, UInt,
	 const std::string, const std::string>(),
	 py::arg("nele"), py::arg("nb_time_steps"),
	 py::arg("dom_size"), py::arg("nu_top"), py::arg("nu_bot"),
	 py::arg("E_top"), py::arg("E_bot"), py::arg("cs_top"),
	 py::arg("cs_bot"), py::arg("tcut_top"), py::arg("tcut_bot"),
	 py::arg("simulation_summary"), py::arg("output_dir")="./",
	 "Constructor for 1D (line) interface between dissimilar bodies")
    .def(py::init< std::vector<UInt>, UInt, std::vector<Real>,
	 Real, Real, Real, UInt,
	 const std::string, const std::string>(),
	 py::arg("nele"), py::arg("nb_time_steps"),
	 py::arg("dom_size"), py::arg("nu"),
	 py::arg("E"), py::arg("cs"), py::arg("tcut"),
	 py::arg("simulation_summary"), py::arg("output_dir")="./",
	 "Constructor for 2D (plane) interface between similar bodies")
    .def(py::init< std::vector<UInt>, UInt, std::vector<Real>,
	 Real, Real, Real, Real, Real, Real, UInt, UInt,
	 const std::string, const std::string>(),
	 py::arg("nele"), py::arg("nb_time_steps"),
	 py::arg("dom_size"), py::arg("nu_top"), py::arg("nu_bot"),
	 py::arg("E_top"), py::arg("E_bot"), py::arg("cs_top"),
	 py::arg("cs_bot"), py::arg("tcut_top"), py::arg("tcut_bot"),
	 py::arg("simulation_summary"), py::arg("output_dir")="./",
	 "Constructor for 2D (plane) interface between dissimilar bodies")
    
    .def("initModel",&SpectralModel::initModel,py::arg("beta")=0,py::arg("blank")=false,
	 "Init the model")
    .def("pauseModel",&SpectralModel::pauseModel,
	 "Pause the model")
    .def("restartModel",&SpectralModel::restartModel,
	 "Restart model from restart files")
    .def("sinusoidalLoading",&SpectralModel::sinusoidalLoading,py::arg("min"),
	 "Set the loading to a sinusoidale shape with the period being the length of the interface")
    .def("readSpatialLoadingFromFile",&SpectralModel::readSpatialLoadingFromFile)
    .def("setLoadingCase",&SpectralModel::setLoadingCase,py::arg("load_in"),py::arg("psi"),py::arg("phi"),py::arg("write")=true)
    .def("setLoadingFromVector",&SpectralModel::setLoadingFromVector)
    .def("setLoadingShape",&SpectralModel::setLoadingShape)
    .def("incrementLoad",&SpectralModel::incrementLoad)
    .def("updateLoads",py::overload_cast<>(&SpectralModel::updateLoads))
    .def("updateLoads",py::overload_cast<Real *>(&SpectralModel::updateLoads))
    .def("initInterfaceFields",&SpectralModel::initInterfaceFields)
    .def("increaseTimeStep",&SpectralModel::increaseTimeStep)	  
    .def("updateDisplacements",&SpectralModel::updateDisplacements)
    .def("setDisplacements",&SpectralModel::setDisplacements)
    .def("setVelocities",&SpectralModel::setVelocities)
    .def("computeInterfaceFields",&SpectralModel::computeInterfaceFields)
    .def("fftOnDisplacements",&SpectralModel::fftOnDisplacements)
    .def("computeStress",&SpectralModel::computeStress)
    .def("printSelfLoad",&SpectralModel::printSelfLoad)
    .def("printSel",&SpectralModel::printSelf)
    // Accessors
    .def("getTime",&SpectralModel::getTime)
    .def("getCurrentTimeStep",&SpectralModel::getCurrentTimeStep)
    .def("getBeta",&SpectralModel::getBeta)
    .def("getDxMin",&SpectralModel::getDxMin)
    .def("getShearWaveSpeeds",&SpectralModel::getShearWaveSpeeds)
    .def("getElementSize",&SpectralModel::getElementSize)
    .def("getNbElements",&SpectralModel::getNbElements)
    .def("getNbTimeSteps",&SpectralModel::getNbTimeSteps)
    .def("getDim",&SpectralModel::getDim)
    .def("getUniformLoading",&SpectralModel::getUniformLoading)
    .def("getInterfaceLaw",&SpectralModel::getInterfaceLaw,py::return_value_policy::reference)
    .def("setInterfaceLaw",&SpectralModel::setInterfaceLaw);
}

} // namespace cracklet
