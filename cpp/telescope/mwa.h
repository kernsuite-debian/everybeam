// mwa.h: Class for MWA telescopes.
// Inherits from Telescope class.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_MWA_H_
#define EVERYBEAM_TELESCOPE_MWA_H_

#include "telescope.h"

#include <casacore/measures/Measures/MPosition.h>

namespace everybeam {

namespace griddedresponse {
class MWAGrid;
class GriddedResponse;
}  // namespace griddedresponse

namespace telescope {

class MWA final : public Telescope {
  friend class griddedresponse::MWAGrid;

 public:
  /**
   * @brief Construct a new MWA object
   *
   * @param ms MeasurementSet
   * @param model Element Response model
   * @param options telescope options
   */
  MWA(const casacore::MeasurementSet &ms, const Options &options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) override;

 private:
  struct MSProperties {
    double delays[16];
    casacore::MPosition array_position;
  };
  MSProperties ms_properties_;
};

}  // namespace telescope
}  // namespace everybeam
#endif  // EVERYBEAM_TELESCOPE_MWA_H_
