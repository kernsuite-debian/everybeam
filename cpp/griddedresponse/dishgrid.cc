// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "dishgrid.h"
#include "../telescope/dish.h"
#include "../circularsymmetric/voltagepattern.h"

#include <aocommon/uvector.h>
#include <aocommon/matrix2x2.h>
#include <aocommon/matrix4x4.h>
#include <aocommon/hmatrix4x4.h>
#include <algorithm>

using aocommon::HMC4x4;
using aocommon::UVector;

namespace everybeam {
namespace griddedresponse {

void DishGrid::Response([[maybe_unused]] BeamMode beam_mode,
                        std::complex<float>* buffer,
                        [[maybe_unused]] double time, double frequency,
                        [[maybe_unused]] size_t station_idx, size_t field_id) {
  const telescope::Dish& dish_telescope =
      static_cast<const telescope::Dish&>(*telescope_);

  double pdir_ra;
  double pdir_dec;
  std::tie(pdir_ra, pdir_dec) = dish_telescope.GetFieldPointing()[field_id];
  const double max_radius_arc_min =
      dish_telescope.GetDishCoefficients()->MaxRadiusInArcMin();
  const double reference_frequency =
      dish_telescope.GetDishCoefficients()->ReferenceFrequency();
  circularsymmetric::VoltagePattern vp(
      dish_telescope.GetDishCoefficients()->GetFrequencies(frequency),
      max_radius_arc_min);
  const aocommon::UVector<double> coefs_vec =
      dish_telescope.GetDishCoefficients()->GetCoefficients(frequency);
  vp.EvaluatePolynomial(coefs_vec, reference_frequency, false);
  vp.Render(buffer, width_, height_, dl_, dm_, ra_, dec_, pdir_ra, pdir_dec,
            l_shift_, m_shift_, frequency);
}

void DishGrid::IntegratedResponse(
    BeamMode /* beam_mode */, float* buffer, double, double frequency,
    size_t field_id, size_t undersampling_factor,
    const std::vector<double>& /*baseline_weights*/) {
  // Copy coordinate members
  const size_t width_original = width_;
  const size_t height_original = height_;
  const double dl_original = dl_;
  const double dm_original = dm_;

  width_ /= undersampling_factor;
  height_ /= undersampling_factor;
  dl_ *= (double(width_original) / double(width_));
  dm_ *= (double(width_original) / double(width_));

  // Init (Hermitian) Mueller matrix for every pixel in the coarse grid
  const size_t npixels = width_ * height_;
  std::vector<HMC4x4> matrices(npixels, HMC4x4::Zero());
  MakeIntegratedDishSnapshot(matrices, frequency, field_id);

  DoFFTResampling(buffer, width_, height_, width_original, height_original,
                  matrices);

  // Reset coordinate members to original values
  width_ = width_original;
  height_ = height_original;
  dl_ = dl_original;
  dm_ = dm_original;
}

void DishGrid::MakeIntegratedDishSnapshot(
    std::vector<aocommon::HMC4x4>& matrices, double frequency,
    size_t field_id) {
  const size_t nstations = telescope_->GetNrStations();
  UVector<std::complex<float>> buffer_undersampled(
      GetStationBufferSize(nstations));
  // Assert that buffer size can accomodate Jones matrix on pixels
  assert(buffer_undersampled.size() >= width_ * height_ * 4);

  ResponseAllStations(BeamMode::kFull, buffer_undersampled.data(), 0.0,
                      frequency, field_id);
  // Loop over the pixels just once, and compute auto correlation
  for (size_t y = 0; y != height_; ++y) {
    for (size_t x = 0; x != width_; ++x) {
      const size_t offset = (y * width_ + x) * 4;
      const aocommon::MC2x2 A(&buffer_undersampled[offset]);
      // Compute Mueller matrix and apply vec trick, see
      // https://en.wikipedia.org/wiki/Kronecker_product#Matrix_equations
      // No need to add (+) and average (*0.5) for auto-correlation
      matrices[y * width_ + x] =
          HMC4x4::KroneckerProduct(A.HermTranspose().Transpose(), A);
    }
  }
}
}  // namespace griddedresponse
}  // namespace everybeam
