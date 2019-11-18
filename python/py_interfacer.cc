#include <interfacer.cc>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cRacklet {
/* -------------------------------------------------------------------------- */

void register_fracture_law_type(py::module&mod){
  py::enum_<FractureLawType>(mod,"FractureLawType")
    .value("_linear_coupled_cohesive",FractureLawType::_linear_coupled_cohesive);
  }
  
//template<FractureLawType F>
void register_interfacer(py::module& mod) {
  FractureLawType const F = _linear_coupled_cohesive;
  py::class_<Interfacer<F>,DataRegister>(mod, "Interfacer")
    .def(py::init<SpectralModel&>())
    .def("createUniformInterface",&Interfacer<F>::createUniformInterface)
    .def("insertPatternfromFile",&Interfacer<F>::insertPatternfromFile)
    .def("createNormalDistributedInterface",&Interfacer<F>::createNormalDistributedInterface)
#ifdef CRACKLET_USE_LIBSURFER
    .def("createBrownianHeterogInterface",&Interfacer<F>::createBrownianHeterogInterface)
#endif
    .def("createThroughArea",&Interfacer<F>::createThroughArea)
    .def("createThroughCrack",&Interfacer<F>::createThroughCrack)
    .def("createThroughWall",&Interfacer<F>::createThroughWall)
    .def("createThroughPolarAsperity",&Interfacer<F>::createThroughPolarAsperity)
    .def("createThroughMultiPolAsperity",&Interfacer<F>::createThroughMultiPolAsperity)
    .def("createThroughCenteredCrack",&Interfacer<F>::createThroughCenteredCrack)
  .def("createThroughLeftSidedCrack",&Interfacer<F>::createThroughLeftSidedCrack)
  .def("createRightPropagatingCrackRoundAsp",&Interfacer<F>::createRightPropagatingCrackRoundAsp)
    .def("createIncohIntfc",&Interfacer<F>::createIncohIntfc);
}

  //void register_interfacer(py::module& mod) {
  
} // namespace cRacklet
