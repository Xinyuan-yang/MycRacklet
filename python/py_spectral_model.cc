#include <spectral_model.hh>
#include <interface_law.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cRacklet {
/* -------------------------------------------------------------------------- */
    
void register_spectral_model(py::module& mod) {      
  py::class_<InterfaceLaw>(mod, "InterfaceLaw");

  py::class_<SpectralModel,DataRegister>(mod, "SpectralModel")
    .def(py::init<>())
    .def(py::init< std::vector<UInt>, UInt, std::vector<Real>,
	 Real, Real, Real, Real, Real, Real, UInt, UInt,
	 const std::string, const std::string>(),
	 py::arg("nele"), py::arg("nb_time_steps"),
	 py::arg("dom_size"), py::arg("nu_top"), py::arg("nu_bot"),
	 py::arg("E_top"), py::arg("E_bot"), py::arg("cs_top"),
	 py::arg("cs_bot"), py::arg("tcut_top"), py::arg("tcut_bot"),
	 py::arg("simulation_summary"), py::arg("output_dir")="./" )
    .def("initModel",&SpectralModel::initModel)
    .def("pauseModel",&SpectralModel::pauseModel)
    .def("restartModel",&SpectralModel::restartModel)
    .def("sinusoidalLoading",&SpectralModel::sinusoidalLoading)

#ifdef CRACKLET_USE_LIBSURFER
    .def("brownianHeterogLoading",&SpectralModel::brownianHeterogLoading)
#endif

    .def("setLoadingCase",&SpectralModel::setLoadingCase)
    .def("updateLoads",py::overload_cast<>(&SpectralModel::updateLoads))
    .def("updateLoads",py::overload_cast<Real *>(&SpectralModel::updateLoads))
    .def("initInterfaceFields",&SpectralModel::initInterfaceFields)
    .def("increaseTimeStep",&SpectralModel::increaseTimeStep)	  
    .def("updateDisplacements",&SpectralModel::updateDisplacements)
    .def("computeInterfaceFields",&SpectralModel::computeInterfaceFields)
    .def("fftOnDisplacements",&SpectralModel::fftOnDisplacements)
    .def("computeStress",&SpectralModel::computeStress)
    .def("printSelfLoad",&SpectralModel::printSelfLoad)
    .def("printSel",&SpectralModel::printSelf)
    // Accessors
    .def("getTime",&SpectralModel::getTime)
    .def("getCurrentTimeStep",&SpectralModel::getCurrentTimeStep)
    .def("getBeta",&SpectralModel::getBeta)
    .def("getElementSize",&SpectralModel::getElementSize)
    .def("getNbElements",&SpectralModel::getNbElements)
    .def("getNbTimeSteps",&SpectralModel::getNbTimeSteps)
    .def("getDim",&SpectralModel::getDim)
    .def("getUniformLoading",&SpectralModel::getUniformLoading)
    .def("getInterfaceLaw",&SpectralModel::getInterfaceLaw,py::return_value_policy::copy);
}

} // namespace cRacklet
