#ifndef __CRACKLET_PY_DATA_REGISTER_HH__
#define __CRACKLET_PY_DATA_REGISTER_HH__

namespace pybind11{
  struct module;
}

namespace cRacklet{

  void register_data_register(pybind11::module & mod);

}

#endif
