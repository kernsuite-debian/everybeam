// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "sphericalharmonicsresponsefixeddirection.h"

#include "common/mathutils.h"
#include "common/sphericalharmonics.h"

namespace everybeam {

SphericalHarmonicsResponseFixedDirection::
    SphericalHarmonicsResponseFixedDirection(
        std::shared_ptr<const SphericalHarmonicsResponse> element_response,
        double theta, double phi)
    : element_response_(std::move(element_response)) {
  assert(element_response_);
  const SphericalHarmonicsResponse::Nms& nms = element_response_->GetNms();

  // Compute all base functions.
  base_functions_.resize(std::array<std::size_t, 2>{nms.shape(0), 2});
  for (size_t i = 0; i < nms.shape(0); ++i) {
    std::tie(base_functions_(i, 0), base_functions_(i, 1)) =
        everybeam::common::F4far_new(nms(i, 2), nms(i, 1), nms(i, 0), theta,
                                     phi);
  }
}

aocommon::MC2x2 SphericalHarmonicsResponseFixedDirection::Response(
    double frequency, double theta, double phi) const {
  if (!element_response_->HasFixedElementIndex()) {
    throw std::runtime_error(
        "SphericalHarmonicsResponseFixedDirection needs an element id, since "
        "SphericalHarmonicsResponse loaded coefficients for all elements.");
  }
  return ComputeResponse(0, frequency);
}

aocommon::MC2x2 SphericalHarmonicsResponseFixedDirection::Response(
    int element_id, double frequency, double theta, double phi) const {
  const std::size_t element_index = element_id;
  if (element_response_->HasFixedElementIndex()) {
    if (element_response_->GetElementIndex() != element_index) {
      throw std::runtime_error("Requested element was not loaded");
    }
  } else {
    if (element_index >= element_response_->GetCoefficients().shape(2)) {
      throw std::runtime_error("Element id is out of range");
    }
  }
  return ComputeResponse(element_index, frequency);
}

aocommon::MC2x2 SphericalHarmonicsResponseFixedDirection::ComputeResponse(
    std::size_t element_index, double frequency) const {
  const SphericalHarmonicsResponse::Coefficients& coefficients =
      element_response_->GetCoefficients();
  const std::size_t frequency_index =
      element_response_->FindFrequencyIndex(frequency);

  aocommon::MC2x2 response = aocommon::MC2x2::Zero();

  for (std::size_t i = 0; i < base_functions_.shape(0); ++i) {
    const std::complex<double> c0 =
        coefficients(0, frequency_index, element_index, i);
    const std::complex<double> c1 =
        coefficients(1, frequency_index, element_index, i);
    const std::complex<double> q2 = base_functions_(i, 0);
    const std::complex<double> q3 = base_functions_(i, 1);
    response[0] += q2 * c0;  // xx
    response[1] += q3 * c0;  // xy
    response[2] += q2 * c1;  // yx
    response[3] += q3 * c1;  // yy
  }

  return response;
}

std::shared_ptr<ElementResponse>
SphericalHarmonicsResponseFixedDirection::FixateDirection(
    const vector3r_t& direction) const {
  const vector2r_t thetaphi = cart2thetaphi(direction);

  return std::make_shared<SphericalHarmonicsResponseFixedDirection>(
      element_response_, thetaphi[0], thetaphi[1]);
}

}  // namespace everybeam