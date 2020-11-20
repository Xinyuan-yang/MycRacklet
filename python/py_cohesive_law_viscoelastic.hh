#ifndef __CRACKLET_PY_COHESIVE_LAW_VISCOELASTIC_HH__
#define __CRACKLET_PY_COHESIVE_LAW_VISCOELASTIC_HH__

namespace pybind11{
  struct module;
}

namespace cRacklet{

  void register_cohesive_law_viscoelastic(pybind11::module & mod);

}

#endif
