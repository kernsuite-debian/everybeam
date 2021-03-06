// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ELEMENTRESPONSE_H
#define EVERYBEAM_ELEMENTRESPONSE_H

#include <complex>
#include <ostream>

#include "common/mutable_ptr.h"

namespace everybeam {

namespace common {
template <typename T>
class MutablePtr;
}

enum ElementResponseModel {
  /// The default will select the default element response model
  /// based on the telescope: e.g. LOFAR will use kHamaker,
  /// OSKAR will select kOSKARSphericalWave.
  kDefault,
  kHamaker,
  kLOBES,
  kOSKARDipole,
  kOSKARSphericalWave
};

std::ostream& operator<<(std::ostream& os, ElementResponseModel model);

/**
 * @brief Virtual class for the element response model. All the
 * (antenna/element) response models inherit from this class.
 *
 */
class ElementResponse {
 public:
  virtual ~ElementResponse() {}

  typedef common::MutablePtr<ElementResponse>
      Ptr;  //!< Pointer to ElementResponse object

  /**
   * @brief Virtual implementation of Response method
   *
   * @param freq Frequency of the plane wave (Hz).
   * @param theta Angle wrt. z-axis (rad)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad)
   * @param result Pointer to 2x2 array of Jones matrix
   */
  virtual void Response(double freq, double theta, double phi,
                        std::complex<double> (&result)[2][2]) const = 0;

  /**
   * @brief Virtual implementation of Response method
   *
   * @param element_id ID of element
   * @param freq Frequency of the plane wave (Hz).
   * @param theta Angle wrt. z-axis (rad)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad)
   * @param result Pointer to 2x2 array of Jones matrix
   */
  virtual void Response(int element_id, double freq, double theta, double phi,
                        std::complex<double> (&result)[2][2]) const {
    Response(freq, theta, phi, result);
  }
};
}  // namespace everybeam
#endif
