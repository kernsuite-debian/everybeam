// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_SPHERICAL_HARMONICS_RESPONSE_H_
#define EVERYBEAM_SPHERICAL_HARMONICS_RESPONSE_H_

#include <cassert>
#include <complex>
#include <optional>
#include <string>
#include <tuple>

#include <xtensor/xtensor.hpp>

#include "elementresponse.h"

namespace everybeam {

/**
 * Computes the response using spherical harmonic coefficients.
 */
class [[gnu::visibility("default")]] SphericalHarmonicsResponse
    : public ElementResponse {
 public:
  using Coefficients = xt::xtensor<std::complex<double>, 4>;
  using Frequencies = xt::xtensor<double, 1>;
  using Nms = xt::xtensor<int, 2>;

  /**
   * Constructor, which loads coefficients for all elements or for a single
   * element.
   * @param coefficients_file An HDF5 file with
   * - coefficients ( 2 x n_frequencies x n_elements x n_nms )
   * - frequencies (n_frequencies)
   * - NMS values (n_nms x 3)
   * @param element_index Index of an element. If specified, only load
   * coefficients for the given index. The index should be less than the number
   * of elements in the HDF5 file.
   * @throw std::runtime_error If the file or the element index is invalid.
   */
  explicit SphericalHarmonicsResponse(
      const std::string& coefficients_file,
      std::optional<std::size_t> element_index = std::nullopt);

  /**
   * Alternate constructor, which allows passing all values using a single
   * argument.
   */
  explicit SphericalHarmonicsResponse(
      std::tuple<std::string, std::optional<std::size_t>> arguments)
      : SphericalHarmonicsResponse(std::get<0>(arguments),
                                   std::get<1>(arguments)) {}

  /**
   * @defgroup Access functions to internal data.
   * @{
   */
  const Coefficients& GetCoefficients() const { return coefficients_; }
  const Frequencies& GetFrequencies() const { return frequencies_; }
  const Nms& GetNms() const { return nms_; }
  bool HasFixedElementIndex() const {
    return static_cast<bool>(element_index_);
  }
  std::size_t GetElementIndex() const {
    assert(element_index_);
    return *element_index_;
  }
  /** @} */

  /** @return The index of the frequency that is closest to 'frequency'. */
  std::size_t FindFrequencyIndex(double frequency) const;

  /**
   * @brief Implementation of the Response overload without element id.
   * @throw std::runtime_error If an element id is required, since this class
   * loaded coefficients for all elements in the constructor.
   */
  aocommon::MC2x2 Response(double frequency, double theta, double phi)
      const override;

  /**
   * @brief Implementation of the Response overload with element id.
   * @throw std::runtime_error If this class only has coefficients for a single
   * element id and the element id does not match.
   */
  aocommon::MC2x2 Response(int element_id, double frequency, double theta,
                           double phi) const override;

  /**
   * Creates an SphericalHarmonicsResponseFixedDirection object with the field
   * quantities (i.e. the basefunctions) for the element response given the
   * direction of interest.
   *
   * @param direction Direction of interest (ITRF, m)
   */
  std::shared_ptr<ElementResponse> FixateDirection(const vector3r_t& direction)
      const final override;

 private:
  /**
   * @brief Does the actual response calculation.
   */
  aocommon::MC2x2 ComputeResponse(std::size_t element_index, double frequency,
                                  double theta, double phi) const;

  Coefficients coefficients_;
  Frequencies frequencies_;
  Nms nms_;

  /// Fixed element index used when constructing this class. If it is valid,
  /// Response() only accepts calls with this index or without index.
  std::optional<std::size_t> element_index_;
};

}  // namespace everybeam

#endif