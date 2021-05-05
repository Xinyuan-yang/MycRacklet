#include <pybind11/pybind11.h>

#ifndef __CRACKLET_PY_RATE_AND_STATE_LAW_HH__
#define __CRACKLET_PY_RATE_AND_STATE_LAW_HH__

namespace cracklet{

  void register_rate_and_state_law(pybind11::module & mod);

}

#endif
