/**
 * @file   py_cracklet.cc
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

namespace cracklet {
  void register_all(pybind11::module & mod) {
    register_interface_law(mod);
    register_data_fields(mod);
    register_data_types(mod);
    register_integrator_types(mod);
    register_data_register(mod);
    register_directions(mod);
    register_spectral_model(mod);
    register_load_control_type(mod);
    register_simulation_driver(mod);
    register_interfacer(mod);
    register_fracture_law_type(mod);
    register_contact_law(mod);
    register_coulomb_law(mod);
    register_reg_coulomb_law(mod);
    register_cohesive_law(mod);
    register_cohesive_law_viscoelastic(mod);
    register_rate_and_state_law(mod);
    
    register_output_format(mod);
    register_data_dumper(mod);
    register_crack_profile(mod);

  }
} // namespace cracklet

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

PYBIND11_MODULE(py11_cracklet, mod){
  mod.doc() = "cracklet python interface";

  cracklet::register_all(mod);
  
} // Module cracklet
