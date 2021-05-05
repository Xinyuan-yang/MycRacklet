#include <pybind11/pybind11.h>

#ifndef __CRACKLET_PY_CONTACT_LAW_HH__
#define __CRACKLET_PY_CONTACT_LAW_HH__

namespace cracklet{

  void register_contact_law(pybind11::module & mod);

}

#endif
