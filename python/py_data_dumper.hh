#include <pybind11/pybind11.h>

#ifndef __CRACKLET_PY_DATA_DUMPER_HH__
#define __CRACKLET_PY_DATA_DUMPER_HH__

namespace cracklet{

  void register_output_format(pybind11::module & mod);
  void register_data_dumper(pybind11::module & mod);

}

#endif
