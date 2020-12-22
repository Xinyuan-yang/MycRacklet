#include "py_interface_law.hh"
#include "py_data_register.hh"
#include "py_spectral_model.hh"
#include "py_simulation_driver.hh"
#include "py_interfacer.hh"
#include "py_contact_law.hh"
#include "py_coulomb_law.hh"
#include "py_regularized_coulomb_law.hh"
#include "py_cohesive_law.hh"
#include "py_cohesive_law_viscoelastic.hh"
#include "py_rate_and_state_law.hh"

#include "py_data_dumper.hh"
#include "py_crack_profile.hh"

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(py_cRacklet, mod){
  mod.doc() = "cRacklet python interface";

  cRacklet::register_interface_law(mod);
  
  cRacklet::register_data_fields(mod);
  cRacklet::register_integrator_types(mod);
  cRacklet::register_data_register(mod);

  cRacklet::register_directions(mod);
  cRacklet::register_spectral_model(mod);
  cRacklet::register_simulation_driver(mod);
  cRacklet::register_interfacer(mod);
  cRacklet::register_fracture_law_type(mod);

  cRacklet::register_contact_law(mod);
  cRacklet::register_coulomb_law(mod);
  cRacklet::register_reg_coulomb_law(mod);
  cRacklet::register_cohesive_law(mod);
  cRacklet::register_cohesive_law_viscoelastic(mod);
  cRacklet::register_rate_and_state_law(mod);

  cRacklet::register_output_format(mod);
  cRacklet::register_data_dumper(mod);
  cRacklet::register_crack_profile(mod);
  
}    
