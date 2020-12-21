#include <pybind11/pybind11.h>

#ifndef __CRACKLET_PY_DATA_REGISTER_HH__
#define __CRACKLET_PY_DATA_REGISTER_HH__

namespace cRacklet{

  void register_data_fields(pybind11::module & mod);
  void register_integrator_types(pybind11::module & mod);
  void register_data_register(pybind11::module & mod);

}

#endif
