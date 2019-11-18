#include <regularized_coulomb_law.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cRacklet {
/* -------------------------------------------------------------------------- */
    
void register_reg_coulomb_law(py::module& mod) {
  py::class_<RegularizedCoulombLaw,ContactLaw>(mod, "RegularizedCoulombLaw")
    .def(py::init<Real,Real,UInt>());
  
}

} // namespace cRacklet
