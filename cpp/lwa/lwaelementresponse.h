// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_LWA_ELEMENTRESPONSE_H_
#define EVERYBEAM_LWA_ELEMENTRESPONSE_H_

#include "../elementresponse.h"

namespace everybeam {

//! Implementation of the LWA response model
class LwaElementResponse : public ElementResponse {
 public:
  LwaElementResponse(){};

  ElementResponseModel GetModel() const override {
    return ElementResponseModel::kLwa;
  }

  aocommon::MC2x2 Response(double freq, double theta,
                           double phi) const override {
    aocommon::MC2x2 response = aocommon::MC2x2::Zero();
    // TODO AST-1385: Add real implementation
    return response;
  };

  // TODO (AST-1385) Add GetInstance() function. Similarly to
  // LOBESElementResponse::GetInstance() this function creates a
  // LwaElementResponse the first time it is invoked, and returns the previously
  // created object in successive invocations. That way, EveryBeam will load the
  // LWA coefficients only once.
};
}  // namespace everybeam
#endif