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
    .def("readSpatialLoadingFromFile",&SpectralModel::readSpatialLoadingFromFile,py::arg("loading_file"),
	 "Read loading from a file")
    .def("setLoadingCase",&SpectralModel::setLoadingCase,py::arg("load_in"),py::arg("psi"),py::arg("phi"),py::arg("write")=true,
	 "Set the loading based on the norm, and the angles phi (in the plane x,y) and psi (in the plane,x,z)")
    .def("setLoadingFromVector",&SpectralModel::setLoadingFromVector,py::arg("loading"),
	 "Set loading from a vector")
    .def("setLoadingShape",&SpectralModel::setLoadingShape,py::arg("shape"),
	 "Spatial enveloppe to modifiy the loading shape")
    .def("incrementLoad",&SpectralModel::incrementLoad,py::arg("increment"),py::arg("loading_direction"),
	 "Increment uniformly the load in a given direction")
    // .def("updateLoads",py::overload_cast<>(&SpectralModel::updateLoads),
    // "Update loading case")
    .def("updateLoads",py::overload_cast<Real *>(&SpectralModel::updateLoads),
	 "Update point-wise loading conidtions using an uniform constant value per dimension")
    .def("initInterfaceFields",&SpectralModel::initInterfaceFields,
	 "Set the initial values of interface fields (strength,traction,velocities) using the interface conditions given in the associated InterfaceLaw")
    .def("increaseTimeStep",&SpectralModel::increaseTimeStep,
	 "Launch the registered computers for current time step and increases time step number")	  
    .def("updateDisplacements",&SpectralModel::updateDisplacements,
	 "Update displacements with velocities")
    .def("setDisplacements",&SpectralModel::setDisplacements,py::arg("displ"),py::arg("side"),
	 "Set displacements fields of one side to a given value")
    .def("setVelocities",&SpectralModel::setVelocities,py::arg("displ"),py::arg("side"),
	 "Set displacements fields of one side to a given value")
    .def("computeInterfaceFields",&SpectralModel::computeInterfaceFields,
	 "Compute interface fields (strength,traction,velocities) using the interface conditions given in the associated InterfaceLaw")
    .def("fftOnDisplacements",&SpectralModel::fftOnDisplacements,
	 "Compute FFT on displacement fields")
    .def("computeStress",&SpectralModel::computeStress,
	 "Compute stresses convolution terms by Backward FFT")
    .def("printSelfLoad",&SpectralModel::printSelfLoad,py::arg("load"),py::arg("psi"),py::arg("phi"),
	 "Dump the load parameters to simulation summary file")
    .def("printSelf",&SpectralModel::printSelf,
	 "Dump the model parameters to simulation summary file")
    // Accessors
    .def("getTime",&SpectralModel::getTime,
	 "Return the current simulation time")
    .def("getCurrentTimeStep",&SpectralModel::getCurrentTimeStep,
	 "Return the current simulation time step")
    .def("getBeta",&SpectralModel::getBeta,
	 "Return stable time step ratio beta")
    .def("getDxMin",&SpectralModel::getDxMin,
	 "Return minimum distance between discretization points")
    .def("getShearWaveSpeeds",&SpectralModel::getShearWaveSpeeds,
	 "Return shear wave speeds of the top and bottom material")
    .def("getElementSize",&SpectralModel::getElementSize,
	 "Return plane discretization")
    .def("getNbElements",&SpectralModel::getNbElements,
	 "Return the number of elements")
    .def("getNbTimeSteps",&SpectralModel::getNbTimeSteps,
	 "Return the number of time steps")
    .def("getDim",&SpectralModel::getDim,
	 "Return model dimension")
    .def("getUniformLoading",&SpectralModel::getUniformLoading,
	 "Return uniform loading vector used to set average interface loading conditions (size=dim)")
    .def("getInterfaceLaw",&SpectralModel::getInterfaceLaw,py::return_value_policy::reference,
	 "Get reference to the FractureLaw")
    .def("setInterfaceLaw",&SpectralModel::setInterfaceLaw,py::arg("itf_law"),
	 "Set reference to the FractureLaw");
}

} // namespace cracklet
