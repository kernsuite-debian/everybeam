// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ELEMENT_RESPONSE_SPHERICAL_WAVE_FIXED_DIRECTION_H_
#define EVERYBEAM_ELEMENT_RESPONSE_SPHERICAL_WAVE_FIXED_DIRECTION_H_

#include <xtensor/xtensor.hpp>

#include "sphericalharmonicsresponse.h"

namespace everybeam {

/**
 * Spherical wave element response model with a fixed direction.
 *
 * Fixing the direction allows reusing the base functions in different
 * Response() calls.
 */
class [[gnu::visibility("default")]] SphericalHarmonicsResponseFixedDirection
    : public ElementResponse {
 public:
  /**
   * @brief Construct a new FixedDirection object that wraps
   *        an existing SphericalHarmonicsResponse object.
   *
   * @param element_response The element response object that should be wrapped.
   * @param theta Fixed theta direction.
   * @param phi Fixed phi direction.
   */
  explicit SphericalHarmonicsResponseFixedDirection(
      std::shared_ptr<const SphericalHarmonicsResponse> element_response,
      double theta, double phi);

  ElementResponseModel GetModel() const final override {
    return element_response_->GetModel();
  }

  /**
   * @brief Override of the Response method without element id.
   * @param frequency Frequency of the plane wave (Hz).
   * @param theta Ignored, since the direction is fixed in the constructor.
   * @param phi Ignored, since the direction is fixed in the constructor.
   * @return The 2x2 Jones matrix.
   */
  aocommon::MC2x2 Response(double frequency, double theta, double phi)
      const final override;

  /**
   * @brief Override of the Response method with element id.
   * @param element_id ID of element
   * @param frequency Frequency of the plane wave (Hz).
   * @param theta Ignored, since the direction is fixed in the constructor.
   * @param phi Ignored, since the direction is fixed in the constructor.
   * @return The 2x2 Jones matrix.
   */
  aocommon::MC2x2 Response(int element_id, double frequency, double theta,
                           double phi) const final override;

  /**
   * Creates an SphericalHarmonicsResponseFixedDirection object with the field
   * quantities (i.e. the basefunctions) for the element response given the
   * direction of interest.
   *
   * This function creates a new object with a fixed direction that can be
   * different from the fixed direction of the current object.
   *
   * @param direction Direction of interest (ITRF, m)
   */
  std::shared_ptr<ElementResponse> FixateDirection(const vector3r_t& direction)
      const final override;

 private:
  aocommon::MC2x2 ComputeResponse(std::size_t element_index, double frequency)
      const;

  std::shared_ptr<const SphericalHarmonicsResponse> element_response_;
  xt::xtensor<std::complex<double>, 2> base_functions_;
};
}  // namespace everybeam

#endif
