// dishgrid.h: Class for computing the circular symmetric (gridded) response.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_DISHGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_DISHGRID_H_

#include "griddedresponse.h"

namespace everybeam {
namespace griddedresponse {

/**
 * @brief Class for computing the gridded response of dish telescopes,
 * e.g. VLA, ATCA.
 *
 */
class DishGrid final : public GriddedResponse {
 public:
  DishGrid(telescope::Telescope* telescope_ptr,
           const coords::CoordinateSystem coordinate_system)
      : GriddedResponse(telescope_ptr, coordinate_system){};

  void CalculateStation(std::complex<float>* buffer, double time,
                        double frequency, size_t station_idx,
                        size_t field_id) override;

  void CalculateAllStations(std::complex<float>* buffer, double time,
                            double frequency, size_t field_id) override;

  virtual void CalculateIntegratedResponse(
      double* buffer, double time, double frequency, size_t field_id,
      size_t undersampling_factor,
      const std::vector<double>& baseline_weights) override;

  virtual void CalculateIntegratedResponse(
      double* buffer, const std::vector<double>& time_array, double frequency,
      size_t field_id, size_t undersampling_factor,
      const std::vector<double>& baseline_weights) override {
    // Time does not play a role in the integrated response of a dish telescope,
    // so call CalculateIntegratedResponse as if it were one time step
    CalculateIntegratedResponse(buffer, 0., frequency, field_id,
                                undersampling_factor, std::vector<double>{0});
  };

 private:
  /**
   * @brief Make integrated snapshot, specialized/simplified for dish
   * telescopes.
   *
   * @param matrices Vector of Mueller matrices
   * @param frequency Frequency (Hz)
   * @param field_id Field id
   */
  void MakeIntegratedDishSnapshot(std::vector<aocommon::HMC4x4>& matrices,
                                  double frequency, size_t field_id);
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_DISHGRID_H_
