#ifndef __CRACKLET_PY_CRACK_PROFILE_HH__
#define __CRACKLET_PY_CRACK_PROFILE_HH__

namespace pybind11{
  struct module;
}

namespace cRacklet{

  void register_crack_profile(pybind11::module & mod);

}

#endif
