#include <contact_law.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cRacklet {
/* -------------------------------------------------------------------------- */
    
void register_contact_law(py::module& mod) {
  py::class_<ContactLaw, std::shared_ptr<ContactLaw>>(mod, "ContactLaw");
      
}

} // namespace cRacklet
