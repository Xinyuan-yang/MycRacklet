#include <pybind11/pybind11.h>

#ifndef __CRACKLET_PY_SIMULATION_DRIVER_HH__
#define __CRACKLET_PY_SIMULATION_DRIVER_HH__

namespace cRacklet{

  void register_load_control_type(pybind11::module & mod);
  void register_simulation_driver(pybind11::module & mod);

}

#endif
