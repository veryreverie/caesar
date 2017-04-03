# ======================================================================
# Provides an include guard for use with CMake.
# ======================================================================

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)

# A macro for include guards in other files.
macro(include_guard name)
  get_property(already_included GLOBAL PROPERTY ${name}_included SET)
  if(${already_included})
    return()
  else()
  endif()
  set_property(GLOBAL PROPERTY ${name}_included true)
endmacro()
