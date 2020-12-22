set(PYCRACKLET_SRCS
  py_cRacklet.cc
  py_interface_law.cc
  py_data_register.cc
  py_spectral_model.cc
  py_simulation_driver.cc
  py_interfacer.cc
  py_contact_law.cc
  py_coulomb_law.cc
  py_regularized_coulomb_law.cc
  py_cohesive_law.cc
  py_cohesive_law_viscoelastic.cc
  py_rate_and_state_law.cc
  py_data_dumper.cc
  py_crack_profile.cc
  )

pybind11_add_module(py_cRacklet ${PYCRACKLET_SRCS})
target_link_libraries(py_cRacklet PUBLIC cRacklet)
set_target_properties(py_cRacklet PROPERTIES DEBUG_POSTIFX "")

set(_python_install_dir
  ${CMAKE_INSTALL_LIBDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)

install(TARGETS py_cRacklet
  LIBRARY DESTINATION ${_python_install_dir})

install(DIRECTORY cRacklet
  DESTINATION ${_python_install_dir}
  FILES_MATCHING PATTERN "*.py")
