#include <cohesive_law.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cRacklet {
/* -------------------------------------------------------------------------- */
    
void register_cohesive_law(py::module& mod) {
  py::class_<CohesiveLaw,InterfaceLaw>(mod, "CohesiveLaw")
    .def(py::init<>())
    .def("initInterfaceConditions",&CohesiveLaw::initInterfaceConditions)
    .def("updateInterfaceConditions",&CohesiveLaw::updateInterfaceConditions)
    .def("restart",&CohesiveLaw::restart)
    .def("preventSurfaceOverlapping",&CohesiveLaw::preventSurfaceOverlapping)

    // Accessor
    
    .def("getNbCohesiveNodes",&CohesiveLaw::getNbCohesiveNodes)
    .def("getNbBrokenNodes",&CohesiveLaw::getNbBrokenNodes);
  
}

} // namespace cRacklet
