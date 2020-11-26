#ifndef __CRACKLET_PY_SPECTRAL_MODEL_HH__
#define __CRACKLET_PY_SPECTRAL_MODEL_HH__

namespace pybind11{
  struct module;
}

namespace cRacklet{

  void register_directions(pybind11::module & mod);
  void register_spectral_model(pybind11::module & mod);

}

#endif
