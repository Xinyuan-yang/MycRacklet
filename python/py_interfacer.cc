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

    case _viscoelastic_coupled_cohesive:
      name = "CohesiveViscoelastic"; 
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
    .value("_viscoelastic_coupled_cohesive",FractureLawType::_viscoelastic_coupled_cohesive)
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
      
      .def("insertPatternfromFile",py::overload_cast<std::string, UInt>(&Interfacer<F>::insertPatternfromFile))
      .def("insertPatternfromFile",py::overload_cast<std::string, std::string,std::string>(&Interfacer<F>::insertPatternfromFile)) // Create an interface based on three files: one for the strength, the other for the opening, one for the residual (If provided)
      .def("createNormalDistributedInterface",&Interfacer<F>::createNormalDistributedInterface)
      .def("createHeterogeneousInterface",&Interfacer<F>::createHeterogeneousInterface,py::arg("crit_nor_opening"),py::arg("max_nor_strength"),py::arg("crit_shr_opening")=std::vector<Real>(),py::arg("max_shr_strength")=std::vector<Real>()) // Create an interface from vectors
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
  
  template<FractureLawType F>
  void wrap_interfacer_viscoelastic(py::module& mod) {
    std::string class_name = FractureLawToString(F);
    py::class_<Interfacer<F>,DataRegister>(mod,class_name.c_str())
      .def(py::init<SpectralModel&>())
      .def("createUniformInterface",&Interfacer<F>::createUniformInterface)
      .def("createThroughArea",&Interfacer<F>::createThroughArea)
      .def("createThroughCrack",&Interfacer<F>::createThroughCrack)
      .def("createThroughWall",&Interfacer<F>::createThroughWall)
      .def("createHeterogeneousInterface",&Interfacer<F>::createHeterogeneousInterface,py::arg("crit_nor_opening"),py::arg("max_nor_strength"),py::arg("crit_shr_opening")=std::vector<Real>(),py::arg("max_shr_strength")=std::vector<Real>()); // Create an interface from vectors
       }
  
  void register_interfacer(py::module& mod) {
    wrap_interfacer_cohesive<_linear_coupled_cohesive>(mod);
    wrap_interfacer_viscoelastic<_viscoelastic_coupled_cohesive>(mod);
    wrap_interfacer_RANDS<_rate_and_state>(mod);
    wrap_interfacer_RANDS<_regularized_rate_and_state>(mod);
    wrap_interfacer_RANDS<_weakening_rate_and_state>(mod);
    
}
  
} // namespace cRacklet
