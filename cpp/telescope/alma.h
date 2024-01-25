// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_ALMA_H_
#define EVERYBEAM_TELESCOPE_ALMA_H_

#include "telescope.h"

#include "../circularsymmetric/coefficients.h"

#include <casacore/measures/Measures/MDirection.h>

namespace everybeam {

namespace griddedresponse {
class DishGrid;
}  // namespace griddedresponse

namespace pointresponse {
class PointResponse;
}  // namespace pointresponse

namespace telescope {

struct AiryParameters {
  AiryParameters(double dish_diameter, double blocked_diameter,
                 double max_radius)
      : dish_diameter_in_m(dish_diameter),
        blocked_diameter_in_m(blocked_diameter),
        maximum_radius_arc_min(max_radius) {}
  double dish_diameter_in_m;
  double blocked_diameter_in_m;
  double maximum_radius_arc_min;
};

/**
 * Provides that ALMA beam pattern, which is implemented as an
 * Airy disk.
 */
class [[gnu::visibility("default")]] Alma final : public Telescope {
 public:
  Alma(const casacore::MeasurementSet& ms, const Options& options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const aocommon::CoordinateSystem& coordinate_system) const override;

  std::unique_ptr<pointresponse::PointResponse> GetPointResponse(double time)
      const override;

  /**
   * @brief Get (ra, dec) pointings of fields.
   *
   * @return Vector of size number of fields, and (ra, dec) pointings as
   * entries.
   */
  const std::vector<std::pair<double, double>>& GetFieldPointing() const {
    return field_pointing_;
  }

  const AiryParameters& GetAiryParameters(size_t station_index) const {
    return parameters_[station_index];
  }

  bool IsHomogeneous() const { return is_homogeneous_; }

 private:
  std::vector<AiryParameters> parameters_;
  /// Store ra, dec pointing per field id from measurement set
  std::vector<std::pair<double, double>> field_pointing_;
  bool is_homogeneous_;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_DISH_H_
