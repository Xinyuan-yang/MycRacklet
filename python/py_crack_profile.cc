#include <crack_profile.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cRacklet {
/* -------------------------------------------------------------------------- */
    
void register_crack_profile(py::module& mod) {
  py::class_<CrackProfile>(mod, "CrackProfile")
    .def(py::init<>())
    .def(py::init<std::vector<UInt>,UInt>());      
}

} // namespace cRacklet
