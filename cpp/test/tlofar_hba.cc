// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "../load.h"
#include "../options.h"
#include "../beammode.h"
#include "../griddedresponse/phasedarraygrid.h"
#include "../pointresponse/phasedarraypoint.h"
#include "../elementresponse.h"
#include "../station.h"
#include "../common/types.h"
#include "../telescope/lofar.h"
#include "../aterms/atermconfig.h"
#include "../aterms/parsetprovider.h"
#include "../msreadutils.h"

#include "config.h"
#include <complex>
#include <cmath>
#include <iostream>
#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>

using aocommon::CoordinateSystem;
using everybeam::ATermSettings;
using everybeam::BeamMode;
using everybeam::BeamNormalisationMode;
using everybeam::ElementResponseModel;
using everybeam::Load;
using everybeam::Options;
using everybeam::ReadTileBeamDirection;
using everybeam::Station;
using everybeam::vector3r_t;
using everybeam::aterms::ATermConfig;
using everybeam::aterms::ParsetProvider;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::PhasedArrayGrid;
using everybeam::pointresponse::PhasedArrayPoint;
using everybeam::pointresponse::PointResponse;
using everybeam::telescope::LOFAR;
using everybeam::telescope::Telescope;

namespace {
// A very simple override of the ParsetProvider. Just falls back on the
// default or on hard-coded values
struct ParsetATerms : public ParsetProvider {
  virtual std::string GetString(
      [[maybe_unused]] const std::string& key) const final override {
    // Not relevant for EveryBeamATerm
    return "";
  }

  std::string GetStringOr([[maybe_unused]] const std::string& key,
                          const std::string& or_value) const final override {
    // Default response model
    return or_value;
  }

  std::vector<std::string> GetStringList(
      [[maybe_unused]] const std::string& key) const final override {
    return std::vector<std::string>{"beam"};
  }

  double GetDoubleOr([[maybe_unused]] const std::string& key,
                     [[maybe_unused]] double or_value) const final override {
    // Update interval set to 1200 (s)
    return 1200.0;
  }
  bool GetBool([[maybe_unused]] const std::string& key) const final override {
    // No use
    return false;
  }
  bool GetBoolOr([[maybe_unused]] const std::string& key,
                 bool or_value) const final override {
    return or_value;
  }
};
}  // namespace

struct HBAFixture {
  HBAFixture() : time(4929192878.008341), frequency(138476562.5) {
    options.element_response_model = ElementResponseModel::kHamaker;
    ms = casacore::MeasurementSet{LOFAR_HBA_MOCK_MS};
    telescope = Load(ms, options);
    coord_system.width = 4;
    coord_system.height = 4;
    coord_system.ra = 2.15374123;
    coord_system.dec = 0.8415521;
    coord_system.dl = 0.5 * M_PI / 180.;
    coord_system.dm = 0.5 * M_PI / 180.;
    coord_system.l_shift = 0.;
    coord_system.m_shift = 0.;
    grid_response = telescope->GetGriddedResponse(coord_system);
    point_response = telescope->GetPointResponse(time);
  }
  ~HBAFixture(){};
  Options options;
  std::unique_ptr<Telescope> telescope;
  std::unique_ptr<GriddedResponse> grid_response;
  std::unique_ptr<PointResponse> point_response;
  aocommon::CoordinateSystem coord_system;

  casacore::MeasurementSet ms;
  double time;
  double frequency;
};

BOOST_FIXTURE_TEST_SUITE(tlofar_hba, HBAFixture)

BOOST_AUTO_TEST_CASE(load_lofar_hba) {
  options.element_response_model = ElementResponseModel::kHamaker;
  // Assert if we indeed have a LOFAR pointer
  BOOST_CHECK(nullptr != dynamic_cast<LOFAR*>(telescope.get()));

  // Assert if correct number of stations
  std::size_t nstations = 70;
  BOOST_CHECK_EQUAL(telescope->GetNrStations(), nstations);

  // Assert if GetStation(stationd_id) behaves properly
  const LOFAR& lofartelescope = static_cast<const LOFAR&>(*telescope.get());
  BOOST_CHECK_EQUAL(lofartelescope.GetStation(0).GetName(), "CS001HBA0");
}

BOOST_AUTO_TEST_CASE(tile_beam_direction) {
  auto lofar_telescope = dynamic_cast<LOFAR*>(telescope.get());
  // Check consistency of the different methods for computing
  // the tile beam direction
  const casacore::MDirection tile_beam_dir_0 =
      lofar_telescope->GetTileBeamDirection();
  const casacore::MDirection tile_beam_dir_1 = ReadTileBeamDirection(ms);
  for (size_t i = 0; i < 3; ++i) {
    BOOST_CHECK_CLOSE(tile_beam_dir_0.getValue().getValue()[i],
                      tile_beam_dir_1.getValue().getValue()[i], 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(element_response) {
  const LOFAR& lofartelescope = static_cast<const LOFAR&>(*telescope.get());
  // Compute and check the Station::ComputeElementResponse for
  // a "randomly selected"  station (station 11)
  // NOTE: this is a regression test in the sense that we only check whether
  // results are consistently reproduced. Reference solution obtained at
  // commit sha 70a286e7dace4616417b0e973a624477f15c9ce3
  //
  // Direction corresponds to one of the itrf directions of the (16) pixels
  // target_element_response is the element response corresponding to this
  // direction
  vector3r_t direction = {0.397408, 0.527527, 0.750855};
  aocommon::MC2x2 target_element_response(
      {-0.164112, -0.000467162}, {-0.843709, -0.00123631},
      {-0.892528, -0.00126278}, {0.0968527, -6.7158e-05});

  const Station& station = lofartelescope.GetStation(11);
  aocommon::MC2x2 element_response =
      station.ComputeElementResponse(time, frequency, direction, false, true);

  // Check whether element_response and target_element_response are "equal"
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_LT(std::abs(element_response[i] - target_element_response[i]),
                   1e-6);
  }

  // Compute station response for station 63 (see also python/test)
  const Station& station63 = lofartelescope.GetStation(63);

  vector3r_t direction_s63 = {0.424588, 0.4629957, 0.7780411};
  vector3r_t station0_dir = {0.4083262, 0.5273447, 0.7451022};
  vector3r_t tile0_dir = {0.4083268, 0.5273442, 0.7451022};
  aocommon::MC2x2 station63_response = station63.Response(
      time, frequency, direction_s63, frequency, station0_dir, tile0_dir);

  aocommon::MC2x2 target_station_response(
      {0.032594235, -0.00023045994}, {0.12204097, -0.00091857865},
      {0.13063535, -0.0010039175}, {-0.029348446, 0.00023882818});

  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_LT(std::abs(station63_response[i] - target_station_response[i]),
                   1e-6);
  }
}

BOOST_AUTO_TEST_CASE(gridded_response) {
  BOOST_CHECK(nullptr != dynamic_cast<PhasedArrayGrid*>(grid_response.get()));
  BOOST_CHECK(nullptr != dynamic_cast<PhasedArrayPoint*>(point_response.get()));

  // Define buffer and get gridded responses
  std::vector<std::complex<float>> antenna_buffer_single(
      grid_response->GetStationBufferSize(1));
  const BeamMode beam_mode = BeamMode::kFull;
  grid_response->Response(beam_mode, antenna_buffer_single.data(), time,
                          frequency, 23, 0);
  BOOST_CHECK_EQUAL(
      antenna_buffer_single.size(),
      std::size_t(coord_system.width * coord_system.height * 2 * 2));
  // LOFARBeam output at pixel (2,2):
  std::vector<std::complex<float>> lofar_p22 = {{-0.175908, -0.000478397},
                                                {-0.845988, -0.00121503},
                                                {-0.89047, -0.00125383},
                                                {0.108123, -5.36076e-05}};

  // Compute response for center pixel via PointResponse
  // One station
  std::complex<float> point_buffer_single_station[4];
  point_response->Response(beam_mode, point_buffer_single_station,
                           coord_system.ra, coord_system.dec, frequency, 23, 0);

  // All stations
  std::vector<std::complex<float>> point_buffer_all_stations(
      point_response->GetAllStationsBufferSize());
  point_response->ResponseAllStations(
      beam_mode, point_buffer_all_stations.data(), coord_system.ra,
      coord_system.dec, frequency, 0);

  // Compare with everybeam
  std::size_t offset_22 = (2 + 2 * coord_system.width) * 4;
  // Offset for station 23
  std::size_t offset_point = 4 * 23;
  for (std::size_t i = 0; i < 4; ++i) {
    // Tolerance is a percentage, so 1e-2 --> 1e-4
    BOOST_CHECK_CLOSE(antenna_buffer_single[offset_22 + i], lofar_p22[i], 1e-2);
    // Following must match exactly, hence, use tighter tolerance
    BOOST_CHECK_CLOSE(point_buffer_single_station[i],
                      antenna_buffer_single[offset_22 + i], 1e-6);
    BOOST_CHECK_CLOSE(point_buffer_all_stations[offset_point + i],
                      antenna_buffer_single[offset_22 + i], 1e-6);
  }

  // LOFARBeam output at pixel (1,3):
  std::vector<std::complex<float>> lofar_p13 = {{-0.158755, -0.000749433},
                                                {-0.816165, -0.00272568},
                                                {-0.863389, -0.00283979},
                                                {0.0936919, 0.000110673}};

  // Compare with everybeam
  std::size_t offset_13 = (1 + 3 * coord_system.width) * 4;
  for (std::size_t i = 0; i < 4; ++i) {
    // Tolerance is a percentage, so 1e-2 --> 1e-4
    BOOST_CHECK_CLOSE(antenna_buffer_single[offset_13 + i], lofar_p13[i], 1e-2);
  }

  // All stations
  std::vector<std::complex<float>> antenna_buffer_all(
      grid_response->GetStationBufferSize(telescope->GetNrStations()));
  grid_response->ResponseAllStations(beam_mode, antenna_buffer_all.data(), time,
                                     frequency, 0);
  BOOST_CHECK_EQUAL(antenna_buffer_all.size(),
                    std::size_t(telescope->GetNrStations() *
                                coord_system.width * coord_system.height * 4));

  // Check consistency of values for station 23
  std::size_t offset_s23 = 23 * coord_system.width * coord_system.height * 4;
  for (std::size_t i = 0; i != antenna_buffer_single.size(); ++i) {
    BOOST_CHECK_LT(
        std::abs(antenna_buffer_all[offset_s23 + i] - antenna_buffer_single[i]),
        1e-6);
  }

  // Check result via aterm calculation
  // Fake the original time for the aterm calculation, accounting for half the
  // update interval
  double time1 = time - 600;
  ATermSettings aterm_settings;
  std::vector<std::complex<float>> aterm_buffer(antenna_buffer_all.size());
  ParsetATerms parset_aterms;
  ATermConfig aterms(telescope->GetNrStations(), coord_system, aterm_settings);
  aterms.Read(ms, parset_aterms);
  aterms.Calculate(aterm_buffer.data(), time1, frequency, 0, nullptr);

  // Check against antenna_buffer_all
  for (std::size_t i = 0; i != aterm_buffer.size(); ++i) {
    BOOST_CHECK_CLOSE(aterm_buffer[i], antenna_buffer_all[i], 1e-6);
  }

  // Save buffer for later reference
  std::vector<std::complex<float>> aterm_ref = aterm_buffer;

  // Result should not change for time increase <1200s
  aterms.Calculate(aterm_buffer.data(), time1 + 1199, frequency, 0, nullptr);
  for (std::size_t i = 0; i != aterm_buffer.size(); ++i) {
    BOOST_CHECK_CLOSE(aterm_buffer[i], aterm_ref[i], 1e-6);
  }

  // Result should change for time increase >=1200s
  aterms.Calculate(aterm_buffer.data(), time1 + 1201, frequency, 0, nullptr);
  for (std::size_t i = 0; i != aterm_buffer.size(); ++i) {
    BOOST_CHECK_GT(std::abs(aterm_buffer[i] - aterm_ref[i]), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(point_response_caching) {
  BOOST_CHECK_EQUAL(point_response->HasTimeUpdate(), true);

  std::vector<std::complex<float>> point_buffer_1(4);
  const BeamMode beam_mode = BeamMode::kFull;
  point_response->Response(beam_mode, point_buffer_1.data(), coord_system.ra,
                           coord_system.dec, frequency, 23, 0);

  BOOST_CHECK_EQUAL(point_response->HasTimeUpdate(), false);

  point_response->UpdateTime(time + 100);
  BOOST_CHECK_EQUAL(point_response->HasTimeUpdate(), true);
  point_response->Response(beam_mode, point_buffer_1.data(), coord_system.ra,
                           coord_system.dec, frequency, 23, 0);

  point_response->SetUpdateInterval(100);
  BOOST_CHECK_EQUAL(point_response->HasTimeUpdate(), true);

  point_response->UpdateTime(time + 100);
  point_response->Response(beam_mode, point_buffer_1.data(), coord_system.ra,
                           coord_system.dec, frequency, 23, 0);

  point_response->UpdateTime(time + 199);
  BOOST_CHECK_EQUAL(point_response->HasTimeUpdate(), false);
  std::vector<std::complex<float>> point_buffer_2(4);
  point_response->Response(beam_mode, point_buffer_2.data(), coord_system.ra,
                           coord_system.dec, frequency, 23, 0);

  for (size_t i = 0; i != point_buffer_1.size(); ++i) {
    BOOST_CHECK_CLOSE(point_buffer_1[i], point_buffer_2[i], 1e-6);
  }

  point_response->UpdateTime(time + 201);
  BOOST_CHECK_EQUAL(point_response->HasTimeUpdate(), true);
  point_response->Response(beam_mode, point_buffer_2.data(), coord_system.ra,
                           coord_system.dec, frequency, 23, 0);

  for (size_t i = 0; i != point_buffer_1.size(); ++i) {
    BOOST_CHECK_PREDICATE(std::not_equal_to<std::complex<float>>(),
                          (point_buffer_1[i])(point_buffer_2[i]));
  }
}

BOOST_AUTO_TEST_CASE(gridded_response_array_factor) {
  // This tests whether "element beam" x "array factor" == "full beam"
  const LOFAR& lofartelescope = static_cast<const LOFAR&>(*telescope.get());
  const Station& station = lofartelescope.GetStation(23);

  // tile0 equals station0
  vector3r_t station0 = {0.408326, 0.527345, 0.745102};
  vector3r_t direction_p13 = {0.397408, 0.527527, 0.750855};
  aocommon::MC2x2 full_beam_p13 = station.Response(
      time, frequency, direction_p13, frequency, station0, station0);
  aocommon::MC2x2 element_beam_p13 = station.ComputeElementResponse(
      time, frequency, direction_p13, false, true);
  aocommon::MC2x2Diag array_factor_p13 = station.ArrayFactor(
      time, frequency, direction_p13, frequency, station0, station0);
  aocommon::MC2x2 full_beam_product_p13 = array_factor_p13 * element_beam_p13;
  for (std::size_t i = 0; i < 4; ++i) {
    // Tolerance is a percentage, so 1e-2 --> 1e-4
    BOOST_CHECK_CLOSE(full_beam_product_p13[i], full_beam_p13[i], 1e-2);
  }

  vector3r_t direction_p22 = {0.408326, 0.527345, 0.745102};
  aocommon::MC2x2 full_beam_p22 = station.Response(
      time, frequency, direction_p22, frequency, station0, station0);
  aocommon::MC2x2 element_beam_p22 = station.ComputeElementResponse(
      time, frequency, direction_p22, false, true);
  aocommon::MC2x2Diag array_factor_p22 = station.ArrayFactor(
      time, frequency, direction_p22, frequency, station0, station0);
  aocommon::MC2x2 full_beam_product_p22 = array_factor_p22 * element_beam_p22;

  for (std::size_t i = 0; i < 4; ++i) {
    // Tolerance is a percentage, so 1e-2 --> 1e-4
    BOOST_CHECK_CLOSE(full_beam_product_p22[i], full_beam_p22[i], 1e-2);
  }
}

BOOST_AUTO_TEST_CASE(differential_beam) {
  // Test with differential beam, single
  Options options_diff_beam = options;
  options_diff_beam.element_response_model = ElementResponseModel::kHamaker;
  options_diff_beam.beam_normalisation_mode = BeamNormalisationMode::kFull;

  // Load (a new) LOFAR Telescope
  std::unique_ptr<Telescope> telescope_diff_beam = Load(ms, options_diff_beam);

  std::unique_ptr<GriddedResponse> grid_response_diff_beam =
      telescope_diff_beam->GetGriddedResponse(coord_system);

  std::vector<std::complex<float>> antenna_buffer_diff_beam(
      grid_response_diff_beam->GetStationBufferSize(1));
  grid_response_diff_beam->Response(
      BeamMode::kFull, antenna_buffer_diff_beam.data(), time, frequency, 15, 0);

  std::size_t offset_22 = (2 + 2 * coord_system.width) * 4;
  double norm_jones_mat = 0.;
  for (std::size_t i = 0; i < 4; ++i) {
    norm_jones_mat += std::norm(antenna_buffer_diff_beam[offset_22 + i]);
  }
  BOOST_CHECK_LT(std::abs(norm_jones_mat - 2.), 1e-6);
}

BOOST_AUTO_TEST_CASE(integrated_beam) {
  // Just check whether IntegratedResponse does run and reproduces
  // results for One time interval
  std::vector<float> antenna_buffer_integrated(
      grid_response->GetIntegratedBufferSize());
  std::vector<double> baseline_weights(
      telescope->GetNrStations() * (telescope->GetNrStations() + 1) / 2, 1.);
  grid_response->IntegratedResponse(BeamMode::kFull,
                                    antenna_buffer_integrated.data(), time,
                                    frequency, 0, 2, baseline_weights);

  // Just check whether some (rather arbitrary) numbers are reproduced
  BOOST_CHECK_LT(std::abs(antenna_buffer_integrated[10] - 0.0309436), 1e-6);
  BOOST_CHECK_LT(std::abs(antenna_buffer_integrated[20] - 0.156267), 1e-6);
  BOOST_CHECK_LT(std::abs(antenna_buffer_integrated.back() - 0.0118085), 1e-6);

  // Two time intervals, should give same output as single time interval
  std::fill(antenna_buffer_integrated.begin(), antenna_buffer_integrated.end(),
            0);
  std::vector<double> tarray = {time, time};
  baseline_weights.resize(baseline_weights.size() * tarray.size());
  std::fill(baseline_weights.begin(), baseline_weights.end(), 1.);
  grid_response->IntegratedResponse(BeamMode::kFull,
                                    antenna_buffer_integrated.data(), tarray,
                                    frequency, 0, 2, baseline_weights);

  BOOST_CHECK_LT(std::abs(antenna_buffer_integrated[10] - 0.0309436), 1e-6);
  BOOST_CHECK_LT(std::abs(antenna_buffer_integrated[20] - 0.156267), 1e-6);
  BOOST_CHECK_LT(std::abs(antenna_buffer_integrated.back() - 0.0118085), 1e-6);

  // Primary beam response on 40 x 40 grid.
  // Validated results were obtained with the following wsclean command
  //
  // wsclean -size 40 40  -scale 900asec -apply-primary-beam  LOFAR_MOCK.ms
  //
  // where LOFAR_MOCK.ms the MS available from
  // https://support.astron.nl/software/ci_data/EveryBeam/L258627-one-timestep.tar.bz2
  //
  // PLEASE NOTE: for the sake of testing, the baseline weights were set to 1
  // in wsclean::lbeamimagemaker, i.e.
  //
  // --- a/lofar/lbeamimagemaker.cpp
  // +++ b/lofar/lbeamimagemaker.cpp
  // @@ -380,7 +380,7 @@ void LBeamImageMaker::makeBeamSnapshot(
  //                                 MC4x4::KroneckerProduct(
  // stationGains[a1].HermTranspose().Transpose(),
  //                                     stationGains[a2]);
  // -          double w = weights.Value(a1, a2);
  // +          double w = 1.;
  //
  std::size_t width_pb = 40, height_pb = 40;
  // (0.25 * M_PI / 180.) equals 900asec
  CoordinateSystem coord_system_pb = {width_pb,
                                      height_pb,
                                      coord_system.ra,
                                      coord_system.dec,
                                      (0.25 * M_PI / 180.),
                                      (0.25 * M_PI / 180.),
                                      coord_system.l_shift,
                                      coord_system.m_shift};

  std::unique_ptr<GriddedResponse> grid_response_pb =
      telescope->GetGriddedResponse(coord_system_pb);

  std::vector<float> antenna_buffer_pb(
      grid_response_pb->GetIntegratedBufferSize());
  std::vector<double> baseline_weights_pb(
      telescope->GetNrStations() * (telescope->GetNrStations() + 1) / 2, 1.);
  grid_response_pb->IntegratedResponse(BeamMode::kFull,
                                       antenna_buffer_pb.data(), time,
                                       frequency, 0, 8, baseline_weights_pb);
  // Check diagonal and off-diagonal term in component 0 and 5 of HMC4x4
  // representation of Mueller matrix
  std::size_t offset_01616 = 16 * width_pb + 16,
              offset_02310 = 23 * width_pb + 10,
              offset_52020 = 5 * width_pb * height_pb + 20 * width_pb + 20,
              offset_51825 = 5 * width_pb * height_pb + 18 * width_pb + 25;

  BOOST_CHECK_LT(std::abs(antenna_buffer_pb[offset_01616] - 0.0203205), 1e-6);
  BOOST_CHECK_LT(std::abs(antenna_buffer_pb[offset_02310] - 0.0111653), 1e-6);
  BOOST_CHECK_LT(std::abs(antenna_buffer_pb[offset_52020] - 0.000154699), 1e-6);
  BOOST_CHECK_LT(std::abs(antenna_buffer_pb[offset_51825] - 0.000133702), 1e-6);
}
BOOST_AUTO_TEST_SUITE_END()
