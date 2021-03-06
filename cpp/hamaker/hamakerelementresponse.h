// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef HAMAKER_ELEMENTRESPONSE_H
#define HAMAKER_ELEMENTRESPONSE_H

#include "../elementresponse.h"
#include "hamakercoeff.h"

#include <memory>

namespace everybeam {

//! Implementation of the Hamaker response model
class HamakerElementResponse : public ElementResponse {
 public:
  virtual void Response(
      double freq, double theta, double phi,
      std::complex<double> (&response)[2][2]) const final override;

  static std::shared_ptr<HamakerElementResponse> GetInstance(
      const std::string &name);

 protected:
  std::string GetPath(const char *) const;

  std::unique_ptr<HamakerCoefficients> coeffs_;
};

class HamakerElementResponseHBA : public HamakerElementResponse {
 public:
  HamakerElementResponseHBA();
};

class HamakerElementResponseLBA : public HamakerElementResponse {
 public:
  HamakerElementResponseLBA();
};

}  // namespace everybeam

#endif
