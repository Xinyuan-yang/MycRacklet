#include <pybind11/pybind11.h>

#ifndef __CRACKLET_PY_COHESIVE_LAW_VISCOELASTIC_HH__
#define __CRACKLET_PY_COHESIVE_LAW_VISCOELASTIC_HH__

namespace cracklet{

  void register_cohesive_law_viscoelastic(pybind11::module & mod);

}

#endif
