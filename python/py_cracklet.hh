/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */

#ifndef PY_CRACKLET_HH_
#define PY_CRACKLET_HH_

namespace cracklet {
void register_all(pybind11::module & mod);
}

#endif /* PY_CRACKLET_HH_ */
