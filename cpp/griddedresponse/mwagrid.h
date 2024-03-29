// mwagrid.h: Class for computing the (gridded) response for the MWA telescope
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_MWAGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_MWAGRID_H_

#include "griddedresponse.h"
#include "../mwabeam/tilebeam2016.h"

#include <memory>

namespace everybeam {
namespace griddedresponse {
class [[gnu::visibility("default")]] MWAGrid final : public GriddedResponse {
 public:
  MWAGrid(const telescope::Telescope* telescope_ptr,
          const aocommon::CoordinateSystem coordinate_system)
      : GriddedResponse(telescope_ptr, coordinate_system){};

  void Response(BeamMode beam_mode, std::complex<float> * buffer, double time,
                double frequency, size_t station_idx, size_t field_id) override;

  void ResponseAllStations(BeamMode beam_mode, std::complex<float> * buffer,
                           double time, double frequency, size_t field_id)
      override {
    HomogeneousAllStationResponse(beam_mode, buffer, time, frequency, field_id);
  }

 private:
  std::unique_ptr<everybeam::mwabeam::TileBeam2016> tile_beam_;

  // Override MakeIntegrated snapshot for efficiency
  void MakeIntegratedSnapshot(BeamMode beam_mode,
                              std::vector<aocommon::HMC4x4> & matrices,
                              double time, double frequency, size_t field_id,
                              const double* baseline_weights_interval) override;
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_MWAGRID_H_
