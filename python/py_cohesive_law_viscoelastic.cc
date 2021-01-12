#include <cohesive_law_viscoelastic.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cRacklet {
/* -------------------------------------------------------------------------- */
    
void register_cohesive_law_viscoelastic(py::module& mod) {
  py::class_<CohesiveLawViscoelastic,CohesiveLaw>(mod, "CohesiveLawViscoelastic")
    .def(py::init<>())
    .def("initLinearFormulation",&CohesiveLawViscoelastic::initLinearFormulation)
    .def("initQuadraticFormulation",&CohesiveLawViscoelastic::initQuadraticFormulation)
    .def("initPowerLawFormulation",&CohesiveLawViscoelastic::initPowerLawFormulation)
    .def("initInterfaceConditions",&CohesiveLawViscoelastic::initInterfaceConditions)
    .def("updateInterfaceConditions",&CohesiveLawViscoelastic::updateInterfaceConditions)
    .def("restart",&CohesiveLawViscoelastic::restart);      
}

} // namespace cRacklet
