#include <crack_profile.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
/* -------------------------------------------------------------------------- */
using namespace pybind11::literals;

namespace py = pybind11;
/* -------------------------------------------------------------------------- */
namespace cracklet {
/* -------------------------------------------------------------------------- */
    
  void register_crack_profile(py::module& mod) {
    py::class_<CrackProfile>(mod, "CrackProfile")
      .def(py::init<>())
      .def(py::init<std::vector<UInt>,UInt>())
      .def("__call__",[](CrackProfile&crack, int i){return crack[i];})
      .def("__call__",[](CrackProfile&crack, py::array_t<UInt> idx){
			py::buffer_info buf = idx.request();
			auto profile = py::array_t<Real>(buf.size);
			py::buffer_info pro = profile.request();
			
			Real *ptr1 = static_cast<Real *>(pro.ptr);
			UInt *ptr2 = static_cast<UInt *>(buf.ptr);
			
			for (size_t i = 0; i < buf.shape[0]; i++){
			  ptr1[i] = crack[ptr2[i]];
			}
			return profile;});      
  }
  
} // namespace cracklet
