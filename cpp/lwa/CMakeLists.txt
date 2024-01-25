# Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(lwa SHARED lwaelementresponse.cc)

install(
  TARGETS lwa
  EXPORT EveryBeamTargets
  DESTINATION lib)

if(DOWNLOAD_LWA)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/coeffs/lwa)
  add_custom_target(
    download_lwa_coefficients
    COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_lwa_coeffs.sh
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/coeffs/lwa)

  # Execute the download script before building the 'lwa' library.
  add_dependencies(lwa download_lwa_coefficients)

  # install coefficients
  message("install LWA coefficients in: " ${EVERYBEAM_DATADIR}/lwa)
  install(FILES "${CMAKE_BINARY_DIR}/coeffs/lwa/LWA_OVRO.h5"
          DESTINATION ${EVERYBEAM_DATADIR})

endif()
