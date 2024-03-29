// griddedresponse.h: Base class for computing the (gridded) telescope response.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_GRIDDEDRESPONSE_H_
#define EVERYBEAM_GRIDDEDRESPONSE_GRIDDEDRESPONSE_H_

#include "./../coords/itrfdirection.h"
#include "./../coords/itrfconverter.h"
#include "../beammode.h"

#include <memory>
#include <vector>
#include <thread>

#include <casacore/measures/Measures/MDirection.h>

#include <aocommon/hmatrix4x4.h>
#include <aocommon/coordinatesystem.h>

namespace everybeam {

namespace telescope {
class Telescope;
}

namespace griddedresponse {

/**
 * @brief Virtual base class to compute the gridded response
 *
 */
class [[gnu::visibility("default")]] GriddedResponse {
 public:
  virtual ~GriddedResponse() {}

  /**
   * @brief Compute the beam for a single station, given a prescribed beam mode.
   * Result is stored in the output buffer, which should allow storing Jones
   * matrices (4 complex floats per pixel).
   *
   * @param beam_mode Selects beam mode (BeamMode::kElement,
   * BeamMode::kArrayFactor or BeamMode::kFull)
   * @param destination Output buffer, compute and set size with
   * GriddedResponse::GetStationBufferSize(1)
   * @param station_idx Station index, must be smaller than number of stations
   * in the Telescope
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param frequency Frequency (Hz)
   */
  virtual void Response(BeamMode beam_mode, std::complex<float> * destination,
                        double time, double freq, size_t station_idx,
                        size_t field_id) = 0;

  /**
   * @brief Compute the array factor for all stations in a Telescope.
   * Result is stored in the output buffer, which should accommodate a Jones
   * matrix (4 complex valued floats) per pixel for each station.
   *
   * @param beam_mode Selects beam mode (element, array factor or full)
   * @param destination Output buffer, compute and set size with
   * GriddedResponse::GetStationBufferSize()
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param frequency Frequency (Hz)
   */
  virtual void ResponseAllStations(
      BeamMode beam_mode, std::complex<float> * destination, double time,
      double frequency, size_t field_id) = 0;

  /**
   * @brief Calculate baseline-integrated (undersampled) beam for a single time
   * step. This function makes use of @ref MakeIntegratedSnapshot(). Subclasses
   * may override MakeIntegratedSnapshot() to implement a more efficient
   * version.
   *
   * @param beam_mode Selects beam mode (element, array factor or full)
   * @param destination Buffer for storing the result, should have size width *
   * height * 16
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param frequency Frequency (Hz)
   * @param field_id Field id
   * @param undersampling_factor Undersampling factor
   * @param baseline_weights Baseline weights, size should equal
   * Telescope::GetNrStations() *  (Telescope::GetNrStations() + 1)/2
   */
  virtual void IntegratedResponse(BeamMode beam_mode, float* destination,
                                  double time, double frequency,
                                  size_t field_id, size_t undersampling_factor,
                                  const std::vector<double>& baseline_weights);

  /**
   * @brief Calculate baseline-integrated beam over multiple time steps.
   * This function makes use of @ref MakeIntegratedSnapshot(). Subclasses
   * may override MakeIntegratedSnapshot() to implement a more efficient
   * version. This function stores all Mueller matrices for the full size
   * image in memory. If this is undesirable, the functions
   * @ref UndersampledIntegratedResponse() and @ref UpsampleResponse() may
   * be used.
   *
   * @param beam_mode Selects beam mode (element, array factor or full)
   * @param destination Buffer for storing the result, should have size width *
   * height * 16
   * @param time_array Vector with probing times, modified Julian date, UTC, in
   * seconds (MJD(UTC), s).
   * @param frequency Frequency (Hz)
   * @param field_id Field id
   * @param undersampling_factor Undersampling factor
   * @param baseline_weights Baseline weights, size should equal
   * (Telescope::GetNrStations() *  (Telescope::GetNrStations() + 1)/2) *
   * time_array.size()
   */
  virtual void IntegratedResponse(
      BeamMode beam_mode, float* destination,
      const std::vector<double>& time_array, double frequency, size_t field_id,
      size_t undersampling_factor, const std::vector<double>& baseline_weights);

  /**
   * Same as @ref IntegratedResponse(), but without performing the upsampling.
   * This function therefore returns the undersampled data.
   * Function @ref UpsampleResponse() can be used to upsample the returned
   * data. This route is useful for big images, for which it is undesirable to
   * hold all 16 elements of the Mueller matrix in memory at the same time.
   * @returns The undersampled data, input to @ref UpsampleResponse().
   */
  virtual std::vector<aocommon::HMC4x4> UndersampledIntegratedResponse(
      BeamMode beam_mode, const std::vector<double>& time_array,
      double frequency, size_t field_id, size_t undersampling_factor,
      const std::vector<double>& baseline_weights);

  /**
   * Upsample a single element from an undersampled response.
   * @param destination Result buffer for one Mueller matrix element of size
   * Width() x Height().
   * @param element_index Value from 0 to 15 indicating which Mueller matrix to
   * upsample.
   * @param undersampled_beam The previously calculated undersampled response
   * @param undersampling_factor Undersampling factor that was used in the
   * response calculation call.
   */
  static void UpsampleResponse(
      float* destination, size_t element_index, size_t width, size_t height,
      const std::vector<aocommon::HMC4x4>& undersampled_beam,
      size_t undersampling_factor);

  std::size_t GetStationBufferSize(std::size_t nstations) const {
    return nstations * width_ * height_ * 4u;
  }

  std::size_t GetIntegratedBufferSize() const { return width_ * height_ * 16u; }

  /**
   * Allow undersampled calculation of the beam grid? If true, the beam
   * will be calculated on a smaller grid and FFT resampled to a larger
   * grid. This is important for arrays that produce large images and for
   * which the beam calculations are expensive, such as MWA, LOFAR and
   * SKA. On the other hand, beams that are fast to calculate and might
   * have sharp cut-offs, should not be undersampled.
   */
  virtual bool PerformUndersampling() const { return true; }

 protected:
  /**
   * @brief Construct a new Gridded Response object
   */
  GriddedResponse(const telescope::Telescope* telescope_ptr,
                  const aocommon::CoordinateSystem& coordinate_system)
      : telescope_(telescope_ptr),
        width_(coordinate_system.width),
        height_(coordinate_system.height),
        ra_(coordinate_system.ra),
        dec_(coordinate_system.dec),
        dl_(coordinate_system.dl),
        dm_(coordinate_system.dm),
        l_shift_(coordinate_system.l_shift),
        m_shift_(coordinate_system.m_shift){};

  static void DoFFTResampling(float* destination, int width_in, int height_in,
                              int width_out, int height_out,
                              const std::vector<aocommon::HMC4x4>& matrices);

  /**
   * Fills the buffer by querying the response of the first antenna and copies
   * it n_stations times. Telescopes with homogeneous antennas can use this
   * function to implement @ref ResponseAllStations().
   *
   * @param buffer should be of size width * height * 4 * n_stations.
   */
  void HomogeneousAllStationResponse(BeamMode beam_mode,
                                     std::complex<float> * buffer, double time,
                                     double frequency, size_t field_id);

  /**
   * Fills the buffer by querying the response of the antennas one by one.
   * Telescopes with heterogeneous antennas can use this function to implement
   * @ref ResponseAllStations(). Homogeneous arrays can use the faster
   * function @ref HomogeneousAllStationResponse() instead.
   *
   * @param buffer should be of size width * height * 4 * n_stations.
   */
  void InhomogeneousAllStationResponse(
      BeamMode beam_mode, std::complex<float> * buffer, double time,
      double frequency, size_t field_id);

  const telescope::Telescope* telescope_;
  size_t width_, height_;
  double ra_, dec_, dl_, dm_, l_shift_, m_shift_;

 private:
  /**
   * @brief Calculate a baseline-integrated snapshot.
   * By default, this function will request the response for all antennas and
   * perform a weighted average. (Partly) homogeneous arrays may implement a
   * faster implementation by overriding this method.
   */
  virtual void MakeIntegratedSnapshot(
      BeamMode beam_mode, std::vector<aocommon::HMC4x4> & matrices, double time,
      double frequency, size_t field_id,
      const double* baseline_weights_interval);
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_GRIDDEDRESPONSE_H_
