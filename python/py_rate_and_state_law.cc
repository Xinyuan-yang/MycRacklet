#include <rate_and_state_law.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cRacklet {
/* -------------------------------------------------------------------------- */
    
void register_rate_and_state_law(py::module& mod) {
  py::class_<RateAndStateLaw,InterfaceLaw>(mod, "RateAndStateLaw")
    .def(py::init<>())
    .def("initStandardFormulation",&RateAndStateLaw::initStandardFormulation)
    .def("initVelocityWeakeningFormulation",&RateAndStateLaw::initVelocityWeakeningFormulation)
    .def("initRegularizedFormulation",&RateAndStateLaw::initRegularizedFormulation)
    .def("setVelocityPredictor",&RateAndStateLaw::setVelocityPredictor)
    .def("initInterfaceConditions",&RateAndStateLaw::initInterfaceConditions)
    .def("updateInterfaceConditions",&RateAndStateLaw::updateInterfaceConditions)
    .def("perturbState",&RateAndStateLaw::perturbState)
    .def("insertPerturbationPatch",&RateAndStateLaw::insertPerturbationPatch)
    .def("insertGaussianPerturbation",&RateAndStateLaw::insertGaussianPerturbation)
    .def("restart",&RateAndStateLaw::restart);
}
} // namespace cRacklet
