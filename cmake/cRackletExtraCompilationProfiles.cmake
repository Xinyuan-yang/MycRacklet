#===============================================================================
# @file   cRackletExtraCompilationProfiles.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Thibault Roch <thibault.roch@epfl.ch>
#
# @date creation: Mon Jan 24 2022
#
# @brief  Compilation profiles
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#===============================================================================

option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." FALSE)
mark_as_advanced(FORCE_COLORED_OUTPUT)
if(FORCE_COLORED_OUTPUT)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_flags(cxx "-fcolor-diagnostics")
  else()
    add_flags(cxx "-fdiagnostics-color=always")
  endif()
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DAKANTU_NDEBUG"
  CACHE STRING "Flags used by the compiler during release builds" FORCE)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG_INIT} -ggdb3"
    CACHE STRING "Flags used by the compiler during debug builds" FORCE)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT} -ggdb3"
    CACHE STRING "Flags used by the compiler during debug builds" FORCE)
endif()

function(declare_compilation_profile name)
  include(CMakeParseArguments)

  cmake_parse_arguments(_args
    "" "COMPILER;LINKER;DOC" "" ${ARGN})

  string(TOUPPER "${name}" _u_name)

  if(NOT _args_DOC)
    string(TOLOWER "${name}" _args_DOC)
  endif()

  if(NOT _args_COMPILER)
    message(FATAL_ERROR "declare_compilation_profile: you should at least give COMPILER flags")
  endif()

  if(NOT _args_LINKER)
    set(_args_LINKER ${_args_COMPILER})
  endif()

  foreach(_flag CXX C Fortran SHARED_LINKER EXE_LINKER)
    set(_stage "compiler")
    set(_flags ${_args_COMPILER})
    if(_stage MATCHES ".*LINKER")
      set(_stage "linker")
      set(_flags ${_args_LINKER})
    endif()
    set(CMAKE_${_flag}_FLAGS_${_u_name} ${_flags}
      CACHE STRING "Flags used by the ${_stage} during coverage builds" FORCE)
    mark_as_advanced(CMAKE_${_flag}_FLAGS_${_u_name})
  endforeach()
endfunction()

# Coverage
declare_compilation_profile(COVERAGE
  COMPILER "-g -ggdb3 -DNDEBUG -DCRACKLET_NDEBUG -O2 --coverage")
