/**
 * @file   py_simulation_driver.hh
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
#include <pybind11/pybind11.h>

#ifndef __CRACKLET_PY_SIMULATION_DRIVER_HH__
#define __CRACKLET_PY_SIMULATION_DRIVER_HH__

namespace cracklet{

  void register_load_control_type(pybind11::module & mod);
  void register_simulation_driver(pybind11::module & mod);

}

#endif
