#include <pybind11/pybind11.h>

#ifndef __CRACKLET_PY_SPECTRAL_MODEL_HH__
#define __CRACKLET_PY_SPECTRAL_MODEL_HH__

namespace cRacklet{

  void register_directions(pybind11::module & mod);
  void register_spectral_model(pybind11::module & mod);

}

#endif
