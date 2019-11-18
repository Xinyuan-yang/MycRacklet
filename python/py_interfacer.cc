#include <interfacer.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cRacklet {
/* -------------------------------------------------------------------------- */
std::string FractureLawToString(FractureLawType F) {

  std::stringstream str;  
  std::string name;
  
  switch(F)
    {
    case _linear_coupled_cohesive:
      name = "LinearCoupledCohesive"; 
      break;
      
    case _rate_and_state:
      name = "RateAndState"; 
      break;

    case _regularized_rate_and_state:
      name = "RegularizedRateAndState"; 
      break;

    case _weakening_rate_and_state:
      name = "WeakeningRateAndState";  
      break;
    }
  
  str << "Interfacer" << name;  
  return str.str();
}
  
void register_fracture_law_type(py::module&mod){
  py::enum_<FractureLawType>(mod,"FractureLawType")
    .value("_linear_coupled_cohesive",FractureLawType::_linear_coupled_cohesive)
    .value("_rate_and_state",FractureLawType::_rate_and_state)
    .value("_regularized_rate_and_state",FractureLawType::_regularized_rate_and_state)
    .value("_weakening_rate_and_state",FractureLawType::_weakening_rate_and_state);
  }
  
template<FractureLawType F>
void wrap_interfacer_RANDS(py::module& mod) {
  std::string class_name = FractureLawToString(F);
  py::class_<Interfacer<F>,DataRegister>(mod,class_name.c_str())
    .def(py::init<SpectralModel&>())
    .def("createUniformInterface",&Interfacer<F>::createUniformInterface);
}

  
template<FractureLawType F>
void wrap_interfacer_cohesive(py::module& mod) {
  std::string class_name = FractureLawToString(F);
  py::class_<Interfacer<F>,DataRegister>(mod,class_name.c_str())
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
    .def("createRightPropagatingCrackRoundAsp",&Interfacer<F>::createRightPropagatingCrackRoundAsp);
    //.def("createIncohIntfc",&Interfacer<F>::createIncohIntfc);
  }
    
  void register_interfacer(py::module& mod) {
    wrap_interfacer_cohesive<_linear_coupled_cohesive>(mod);
    wrap_interfacer_RANDS<_rate_and_state>(mod);
    wrap_interfacer_RANDS<_regularized_rate_and_state>(mod);
    wrap_interfacer_RANDS<_weakening_rate_and_state>(mod);

}
  
} // namespace cRacklet
