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
    
void register_interface_law(py::module& mod) {
  py::class_<InterfaceLaw>(mod, "InterfaceLaw");
      
}

} // namespace cracklet
