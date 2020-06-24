#ifndef __CRACKLET_PY_INTERFACE_LAW_HH__
#define __CRACKLET_PY_INTERFACE_LAW_HH__

namespace pybind11{
  struct module;
}

namespace cRacklet{

  void register_interface_law(pybind11::module & mod);

}

#endif
