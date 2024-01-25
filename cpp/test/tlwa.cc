// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../telescope/lwa.h"

#include "../load.h"
#include <boost/test/unit_test.hpp>
#include "config.h"
#include "../correctionmode.h"

BOOST_AUTO_TEST_SUITE(lwa)

BOOST_AUTO_TEST_CASE(get_telescope_type) {
  casacore::MeasurementSet ms(LWA_MOCK_PATH);
  BOOST_ASSERT(everybeam::GetTelescopeType(ms) ==
               everybeam::TelescopeType::kOvroLwaTelescope);
}

BOOST_AUTO_TEST_CASE(load) {
  // Load LWA Telescope
  everybeam::Options options;
  options.element_response_model = everybeam::ElementResponseModel::kLwa;
  std::unique_ptr<everybeam::telescope::Telescope> telescope =
      everybeam::Load(LWA_MOCK_PATH, options);

  // Check that the number of station is correct.
  BOOST_CHECK_EQUAL(telescope->GetNrStations(), size_t{352});

  // Check that we have an LWA pointer.
  const everybeam::telescope::Lwa* lwa =
      dynamic_cast<const everybeam::telescope::Lwa*>(telescope.get());
  BOOST_REQUIRE(lwa);

  // Check that the first 5 antennas have the correct name as in the MS.
  BOOST_CHECK_EQUAL(lwa->GetStation(0).GetName(), "LWA266");
  BOOST_CHECK_EQUAL(lwa->GetStation(1).GetName(), "LWA259");
  BOOST_CHECK_EQUAL(lwa->GetStation(2).GetName(), "LWA268");
  BOOST_CHECK_EQUAL(lwa->GetStation(3).GetName(), "LWA267");
  BOOST_CHECK_EQUAL(lwa->GetStation(4).GetName(), "LWA271");

  // Check that frequency information is read correctly.
  BOOST_CHECK_EQUAL(lwa->GetNrChannels(), 1u);
  BOOST_CHECK_EQUAL(lwa->GetChannelFrequency(0), 59335937.5);
}

BOOST_AUTO_TEST_SUITE_END()
