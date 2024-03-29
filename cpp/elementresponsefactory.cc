// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// elementresponse.cc describes why ElementResponse::GetInstance is in this
// separate file.

#include "elementresponse.h"

#include <iostream>

#include "hamaker/hamakerelementresponse.h"
#include "oskar/oskarelementresponse.h"
#include "lobes/lobeselementresponse.h"
#include "lwa/lwaelementresponse.h"
#include "options.h"

namespace everybeam {

std::shared_ptr<const ElementResponse> ElementResponse::GetInstance(
    ElementResponseModel model, const std::string& name,
    const Options& options) {
  switch (model) {
    case ElementResponseModel::kHamaker:
      return std::make_shared<HamakerElementResponse>(name);
    case ElementResponseModel::kHamakerLba:
      return std::make_shared<HamakerElementResponse>("LBA");
    case ElementResponseModel::kOSKARDipole:
      return std::make_shared<OSKARElementResponseDipole>();
    case ElementResponseModel::kOSKARSphericalWave:
      return std::make_shared<OSKARElementResponseSphericalWave>();
    case ElementResponseModel::kLOBES:
      try {
        return LOBESElementResponse::GetInstance(name, options);
      } catch (const std::runtime_error& e) {
        static bool first_time = true;
        if (first_time) {
          std::cout << "Creating LOBESElementResponse for station " << name
                    << " failed because: \n"
                    << e.what() << '\n'
                    << "Switching to HamakerElementResponse instead. Further "
                       "warnings of this type will be suppressed.\n";
          first_time = false;
        }
        return std::make_shared<HamakerElementResponse>(name);
      }
    case ElementResponseModel::kLwa:
      // TODO AST-1385: use LwaElementResponse::GetInstance() instead
      return std::make_shared<LwaElementResponse>();
    default:
      std::stringstream message;
      message << "The requested element response model '" << model
              << "' is not implemented.";
      throw std::runtime_error(message.str());
  }
}

}  // namespace everybeam
