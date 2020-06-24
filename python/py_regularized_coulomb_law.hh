#ifndef __CRACKLET_PY_REG_COULOMB_LAW_HH__
#define __CRACKLET_PY_REG_COULOMB_LAW_HH__

namespace pybind11{
  struct module;
}

namespace cRacklet{

  void register_reg_coulomb_law(pybind11::module & mod);

}

#endif
