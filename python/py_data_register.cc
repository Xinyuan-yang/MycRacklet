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
    
void register_data_register(py::module& mod) {
  py::class_<DataRegister>(mod, "DataRegister")
    .def(py::init<>())
    .def("registerParameter",&DataRegister::registerParameter)
    .def("readInputFile",&DataRegister::readInputFile)
    .def("getParameter",&DataRegister::getParameter)
    .def("registerComputer",&DataRegister::registerComputer)
    .def("getComputer",&DataRegister::getComputer)
    .def("getCrackTipPosition",&DataRegister::getCrackTipPosition);

}

} // namespace cRacklet
