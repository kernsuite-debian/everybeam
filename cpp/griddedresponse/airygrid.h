// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_AIRYGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_AIRYGRID_H_

#include "griddedresponse.h"

namespace everybeam {
namespace griddedresponse {

/**
 * @brief Class for computing the gridded response of (possibly heterogenous)
 * telescopes with an Airy disk response, e.g. ALMA.
 */
class [[gnu::visibility("default")]] AiryGrid final : public GriddedResponse {
 public:
  AiryGrid(const telescope::Telescope* telescope_ptr,
           const aocommon::CoordinateSystem coordinate_system)
      : GriddedResponse(telescope_ptr, coordinate_system){};

  void Response(BeamMode beam_mode, std::complex<float> * buffer, double time,
                double frequency, size_t station_idx, size_t field_id) override;

  void ResponseAllStations(BeamMode beam_mode, std::complex<float> * buffer,
                           double time, double frequency, size_t field_id)
      override;

  void IntegratedResponse(
      BeamMode beam_mode, float* buffer,
      [[maybe_unused]] const std::vector<double>& time_array, double frequency,
      size_t field_id, size_t undersampling_factor,
      [[maybe_unused]] const std::vector<double>& baseline_weights) override {
    // The integrated response of an Airy disk telescope is independent of time,
    // so call IntegratedResponse as if it were one time step
    GriddedResponse::IntegratedResponse(beam_mode, buffer, 0.0, frequency,
                                        field_id, undersampling_factor, {0.0});
  }

  bool PerformUndersampling() const override { return false; }

 private:
  /**
   * @brief Make integrated snapshot, specialized/simplified for dish
   * telescopes.
   *
   * @param matrices Vector of Mueller matrices
   * @param frequency Frequency (Hz)
   * @param field_id Field id
   */
  void MakeIntegratedDishSnapshot(std::vector<aocommon::HMC4x4> & matrices,
                                  double frequency, size_t field_id);
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_DISHGRID_H_
