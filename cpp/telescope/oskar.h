// oskar.h: Base class for computing the response for the OSKAR
// telescope.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_OSKAR_H_
#define EVERYBEAM_TELESCOPE_OSKAR_H_

#include "../station.h"
#include "../elementresponse.h"
#include "phasedarray.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <memory>

namespace everybeam {

namespace griddedresponse {
class OSKARGrid;
class GriddedResponse;
}  // namespace griddedresponse

namespace telescope {

//! OSKAR telescope class
class OSKAR final : public PhasedArray {
  friend class griddedresponse::OSKARGrid;

 public:
  /**
   * @brief Construct a new OSKAR object
   *
   * @param ms MeasurementSet
   * @param model Element Response model
   * @param options telescope options
   */
  OSKAR(const casacore::MeasurementSet &ms, const Options &options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) override;

  //! Get the tile beam direction, equal to delay direction for OSKAR!
  virtual casacore::MDirection GetTileBeamDirection() const final override {
    std::cout << "OSKAR has no tile. tile_beam_dir is equal to the delay_dir."
              << std::endl;
    return ms_properties_.tile_beam_dir;
  };

  //! Get the preapplied beam direction, equal to delay direction for OSKAR!
  virtual casacore::MDirection GetPreappliedBeamDirection()
      const final override {
    std::cout << "OSKAR has no preapplied beam direction (yet). "
                 "preapplied_beam_dir is equal to the delay_dir."
              << std::endl;
    return ms_properties_.preapplied_beam_dir;
  };
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_OSKAR_H_
