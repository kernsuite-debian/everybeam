#include "alma.h"
#include "../griddedresponse/airygrid.h"
#include "../pointresponse/airypoint.h"

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using everybeam::griddedresponse::AiryGrid;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::pointresponse::AiryPoint;
using everybeam::pointresponse::PointResponse;

namespace everybeam::telescope {

Alma::Alma(const casacore::MeasurementSet& ms,
           const everybeam::Options& options)
    : Telescope(ms, options) {
  casacore::MSField field_table = ms.field();
  casacore::ArrayColumn<double> pointing_dir_col(
      field_table, casacore::MSField::columnName(casacore::MSField::DELAY_DIR));

  for (std::size_t field_id = 0; field_id != field_table.nrow(); ++field_id) {
    casacore::Array<double> pdir = pointing_dir_col(field_id);
    double pdir_ra = *pdir.cbegin();
    double pdir_dec = *(pdir.cbegin() + 1);
    field_pointing_.emplace_back(pdir_ra, pdir_dec);
  }

  casacore::MSAntenna antenna_table = ms.antenna();
  casacore::ScalarColumn<double> diameter_col(
      antenna_table,
      casacore::MSAntenna::columnName(casacore::MSAntenna::DISH_DIAMETER));
  const size_t n_antennas = diameter_col.nrow();
  parameters_.reserve(n_antennas);
  is_homogeneous_ = true;
  int first_diameter = std::round(diameter_col(0));
  for (size_t i = 0; i != n_antennas; ++i) {
    int rounded_diameter = std::round(diameter_col(i));
    is_homogeneous_ = is_homogeneous_ && (first_diameter == rounded_diameter);
    /*
     * ALMA's parameters: (dish size, blocked size, maximum distance to compute)
     * 12M:  10.7, 0.75, 1.784
     * 7M:  6.25, 0.75, 3.568
     * These values were provided by Dirk Petry and Tony Mroczkowski, and can
     * also be found in Casa (source file
     * synthesis/TransformMachines/PBMath.cc). The maximum distance is provided
     * in degrees at GHz, whereas the Airy disk code requires arcmin at GHz,
     * hence the factor of 60.
     */
    if (rounded_diameter == 12) {
      parameters_.emplace_back(10.7, 0.75, 1.784 * 60.0);
    } else if (rounded_diameter == 7) {
      parameters_.emplace_back(6.25, 0.75, 3.568 * 60.0);
    } else
      throw std::runtime_error(
          "Diameter of ALMA dish is set to neither 12 m or 7 m: only these two "
          "diameters are supported");
  }

  SetIsTimeRelevant(false);
}

std::unique_ptr<GriddedResponse> Alma::GetGriddedResponse(
    const aocommon::CoordinateSystem& coordinate_system) const {
  return std::make_unique<AiryGrid>(this, coordinate_system);
}

std::unique_ptr<PointResponse> Alma::GetPointResponse(double time) const {
  return std::make_unique<AiryPoint>(this, time);
}

}  // namespace everybeam::telescope
