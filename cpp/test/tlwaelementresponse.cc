// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../lwa/lwaelementresponse.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(lwaelementresponse)

BOOST_AUTO_TEST_CASE(get_model) {
  auto lwa_reponse = std::make_shared<everybeam::LwaElementResponse>();
  BOOST_TEST(lwa_reponse->GetModel() == everybeam::ElementResponseModel::kLwa);
}

BOOST_AUTO_TEST_SUITE_END()
