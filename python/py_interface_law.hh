#include <pybind11/pybind11.h>

#ifndef __CRACKLET_PY_INTERFACE_LAW_HH__
#define __CRACKLET_PY_INTERFACE_LAW_HH__

namespace cRacklet{

  void register_interface_law(pybind11::module & mod);

}

#endif
