/**
 * @file   py_interfacer.cc
 * @author Thibault Roch <thibault.roch@epfl.ch>
 *
 * @section LICENSE
 *
 * cRacklet - A spectral boundary integral method for interface fracture simulation
 * Copyright (©) 2012 - 2013 Fabian Barras
 *               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 *               Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 * 
 * cRacklet is free software: you can redistribute it and/or modify it under the terms 
 * of the GNU General Public License as published by the Free Software Foundation, 
 * either version 3 of the License, or (at your option) any later version.
 * 
 * cRacklet is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program.  
 * If not, see <http://www.gnu.org/licenses/>.
 */
/* -------------------------------------------------------------------------- */
#include <interfacer.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cracklet {
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

    case _regularized_weakening_rate_and_state:
      name = "RegularizedWeakeningRateAndState"; 
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
    .value("_regularized_weakening_rate_and_state",FractureLawType::_regularized_weakening_rate_and_state)
    .value("_weakening_rate_and_state",FractureLawType::_weakening_rate_and_state);
  }
  
  template<FractureLawType F>
  void wrap_interfacer_RANDS(py::module& mod) {
    std::string class_name = FractureLawToString(F);
    py::class_<Interfacer<F>,DataRegister>(mod,class_name.c_str())
      .def(py::init<SpectralModel&>())
      .def("createUniformInterface",&Interfacer<F>::createUniformInterface,"Create a uniform layer on the entire interface. Required parameters should be registered through the simulation parameters.");
  }
  
  
  template<FractureLawType F>
  void wrap_interfacer_cohesive(py::module& mod) {
    std::string class_name = FractureLawToString(F);
    py::class_<Interfacer<F>,DataRegister>(mod,class_name.c_str())
      .def(py::init<SpectralModel&>())
      .def("createUniformInterface",&Interfacer<F>::createUniformInterface,"Create a uniform layer on the entire interface. Required parameters should be registered through the simulation parameters.")
      
      .def("insertPatternfromFile",[](Interfacer<F> & self, std::string filename, UInt origin=0) {self.insertPatternfromFile(filename,origin);
				   },
	"Tune interface properties from a text file starting at a given position")
      .def("insertPatternfromFile",[](Interfacer<F> & self, std::string file_strength, std::string file_opening, std::string file_residual = "None") {self.insertPatternfromFile(file_strength,file_opening,file_residual);
				   },
	"Create an interface from two text files, one for the critical stress the other for the critical opening, and optionally one for the residual strength")
      .def("createNormalDistributedInterface",&Interfacer<F>::createNormalDistributedInterface,
	   py::arg("crit_nor_opening"),py::arg("crit_shr_opening"),py::arg("max_nor_strength"),py::arg("max_shr_strength"),py::arg("stddev"),py::arg("seed"), 
	   "Create an heterogeneous interface following normal distribution of strength")
      .def("createHeterogeneousInterface",&Interfacer<F>::createHeterogeneousInterface,py::arg("crit_nor_opening"),py::arg("max_nor_strength"),py::arg("res_nor_strength")=std::vector<Real>(),py::arg("crit_shr_opening")=std::vector<Real>(),py::arg("max_shr_strength")=std::vector<Real>(),py::arg("res_shr_strength")=std::vector<Real>(), "Create an heterogeneous interface from vectors") // Create an interface from vectors
      .def("createThroughArea",&Interfacer<F>::createThroughArea,
	   py::arg("area_start"),py::arg("area_end"),py::arg("cracking_index"),py::arg("ratio_max_nor_strength")=1,py::arg("ratio_max_shr_strength")=1,py::arg("ratio_max_crit_nor_opening")=1,py::arg("ratio_max_crit_shr_opening")=1,py::arg("variation_rather_than_ratio")=0, 
	   "create a z-invariant(=through) area between x=start and x=end of given cracking_index and with properties given by new_prop = ratio*current_prop, if variation_rather_than_ratio=0, new_prop = ratio+current_prop, if variation_rather_than_ratio=1")
      .def("createThroughCrack",&Interfacer<F>::createThroughCrack,
	   py::arg("crack_start"),py::arg("crack_end"),
	   "Create a crack between x=crack_start and x=crack_end ")
      .def("createThroughWall",&Interfacer<F>::createThroughWall,
	   py::arg("wall_start"),py::arg("wall_end"),
	   "Create a wall (large thoughness) between x=wall_start and x=wall_end")
      .def("createThroughPolarAsperity",&Interfacer<F>::createThroughPolarAsperity,
	   py::arg("position"),py::arg("width"),py::arg("delta_max_nor_strength"),py::arg("delta_max_shr_strength"),py::arg("delta_crit_nor_opening"),py::arg("delta_crit_shr_opening"),py::arg("polarity"),
	   "Create an asperity made of weaker(-delta) then stronger(+delta) areas at given position and given width")
      .def("createThroughMultiPolAsperity",&Interfacer<F>::createThroughMultiPolAsperity,
	   py::arg("start"),py::arg("end"),py::arg("number"),py::arg("delta_max_nor_strength"),py::arg("delta_max_shr_strength"),py::arg("delta_crit_nor_opening"),py::arg("delta_crit_shr_opening"),py::arg("polarity"),
	   "Create an interface made of multiple weaker(-delta) and stronger(+delta) areas"
	   )
      .def("createThroughCenteredCrack",&Interfacer<F>::createThroughCenteredCrack,
	   py::arg("initial_crack_size"),py::arg("crit_nor_opening"),py::arg("crit_shr_opening"),py::arg("max_nor_strength"),py::arg("max_shr_strength"),
	   "create an interface with a centered crack")
      .def("createThroughLeftSidedCrack",&Interfacer<F>::createThroughLeftSidedCrack,
	   py::arg("initial_crack_size"),py::arg("crit_nor_opening"),py::arg("crit_shr_opening"),py::arg("max_nor_strength"),py::arg("max_shr_strength"),
	   "Create an interface with a left-sided crack")
      .def("createRightPropagatingCrackRoundAsp",&Interfacer<F>::createRightPropagatingCrackRoundAsp,
	   py::arg("initial_crack_size"),py::arg("crit_nor_opening"),py::arg("crit_shr_opening"),py::arg("max_nor_strength"),py::arg("max_shr_strength"),py::arg("radius"),py::arg("asp_ctr"),py::arg("ratio_strength"),py::arg("ratio_crit_open")=1,   
	   "create right propagating through crack meeting a circular asperity whose strength = ratio_strength*interface_strength and crit_opening = ratio_crit_open*interface_crit_opening"
	   );
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
      .def("createHeterogeneousInterface",&Interfacer<F>::createHeterogeneousInterface,py::arg("crit_nor_opening"),py::arg("max_nor_strength"),py::arg("res_nor_strength")=std::vector<Real>(),py::arg("crit_shr_opening")=std::vector<Real>(),py::arg("max_shr_strength")=std::vector<Real>(),py::arg("res_shr_strength")=std::vector<Real>()); // Create an interface from vectors
      }
  
  void register_interfacer(py::module& mod) {
    wrap_interfacer_cohesive<_linear_coupled_cohesive>(mod);
    wrap_interfacer_viscoelastic<_viscoelastic_coupled_cohesive>(mod);
    wrap_interfacer_RANDS<_rate_and_state>(mod);
    wrap_interfacer_RANDS<_regularized_rate_and_state>(mod);
    wrap_interfacer_RANDS<_regularized_weakening_rate_and_state>(mod);
    wrap_interfacer_RANDS<_weakening_rate_and_state>(mod);
    
}
  
} // namespace cracklet
