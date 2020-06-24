#ifndef __CRACKLET_PY_CONTACT_LAW_HH__
#define __CRACKLET_PY_CONTACT_LAW_HH__

namespace pybind11{
  struct module;
}

namespace cRacklet{

  void register_contact_law(pybind11::module & mod);

}

#endif
