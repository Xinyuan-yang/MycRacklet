#ifndef __CRACKLET_PY_INTERFACER_HH__
#define __CRACKLET_PY_INTERFACER_HH__

#include <iostream>

namespace pybind11{
  struct module;
}

namespace cRacklet{

  void register_interfacer(pybind11::module & mod);
  void register_fracture_law_type(pybind11::module & mod);
  //std::string FractureLawToString(FractureLawType F);
  
}

#endif
