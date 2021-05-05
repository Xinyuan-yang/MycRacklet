#include <coulomb_law.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cracklet {
/* -------------------------------------------------------------------------- */
    
void register_coulomb_law(py::module& mod) {
  py::class_<CoulombLaw,ContactLaw,std::shared_ptr<CoulombLaw>>(mod, "CoulombLaw")
    .def(py::init<Real>());
  
}

} // namespace cracklet
