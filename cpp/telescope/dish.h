// dish.h: Base class for dish telescopes (VLA, ATCA, ...).
// Inherits from Telescope class.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_DISH_H_
#define EVERYBEAM_TELESCOPE_DISH_H_

#include "telescope.h"

#include <casacore/measures/Measures/MDirection.h>

namespace everybeam {

namespace griddedresponse {
class DishGrid;
class GriddedResponse;
}  // namespace griddedresponse

namespace telescope {

/**
 * This class calculates the a-terms for dishes with a circularly symmetric
 * response.
 */
class Dish : public Telescope {
  friend class griddedresponse::DishGrid;

 public:
  Dish(const casacore::MeasurementSet &ms, const Options &options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) override;

 private:
  struct MSProperties {
    std::vector<std::pair<double, double>> field_pointing;
  };
  MSProperties ms_properties_;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_DISH_H_
