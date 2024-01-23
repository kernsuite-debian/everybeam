# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(
  aterms OBJECT
  atermconfig.cc
  atermbase.cc
  atermresampler.cc
  everybeamaterm.cc
  fitsatermbase.cc
  fitsaterm.cc
  dldmaterm.cc
  pafbeamterm.cc
  h5parmaterm.cc
  klfitter.cc
  klfittingaterm.cc
  fourierfitter.cc
  fourierfittingaterm.cc
  cache.cc)

# The aterms object library is included in the SHARED library everybeam
# Compilation must therefore generate relocatable code
set_property(TARGET aterms PROPERTY POSITION_INDEPENDENT_CODE ON)

# Link to XTensor libraries, which are loaded using FetchXTensor.cmake.
target_link_libraries(aterms xtensor xtl xtensor-blas xtensor-fftw)

install(
  FILES atermconfig.h
        parsetprovider.h
        # May not need to expose these
        cache.h
        atermbase.h
        atermbeam.h
        atermresampler.h
        #
        fitsatermbase.h
        fitsaterm.h
        everybeamaterm.h
        dldmaterm.h
        pafbeamterm.h
        h5parmaterm.h
        klfittingaterm.h
        fourierfittingaterm.h
  DESTINATION "include/${CMAKE_PROJECT_NAME}/aterms")
