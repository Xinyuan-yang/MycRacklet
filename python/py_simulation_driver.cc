#include <simulation_driver.cc>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cRacklet {
  /* -------------------------------------------------------------------------- */
  
  void register_load_control_type(py::module& mod) {
    py::enum_<LoadControlType>(mod,"LoadControlType")
      .value("_space_control",LoadControlType::_space_control)
      .value("_time_control",LoadControlType::_time_control)
      .export_values();
  }
  
  void register_simulation_driver(py::module& mod) {
    py::class_<SimulationDriver,DataRegister>(mod, "SimulationDriver")
      .def(py::init<SpectralModel&, Real>(),py::arg("model_to_drive"),py::arg("beta")=0.0)
      .def(py::init<SpectralModel&, Real, Real>())
      .def("initConstantLoading",&SimulationDriver::initConstantLoading)
      .def("initLoadingFromFile",&SimulationDriver::initLoadingFromFile)
      .def("initConstantSpeed",&SimulationDriver::initConstantSpeed)
      .def("solveStep",&SimulationDriver::solveStep)
      .def("writeLoading",&SimulationDriver::writeLoading)
      .def("launchCrack",&SimulationDriver::launchCrack,py::arg("crack_start"),py::arg("launched_size"),py::arg("v_init"),py::arg("one_side_propagation")=true);
  }
  
} // namespace cRacklet
