// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "load.h"

#include "telescope/alma.h"
#include "telescope/dish.h"
#include "telescope/lofar.h"
#include "telescope/lwa.h"
#include "telescope/mwa.h"
#include "telescope/oskar.h"
#include "telescope/skamid.h"

#include "circularsymmetric/atcacoefficients.h"
#include "circularsymmetric/gmrtcoefficients.h"
#include "circularsymmetric/meerkatcoefficients.h"
#include "circularsymmetric/vlacoefficients.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <string>
#include <stdexcept>

namespace everybeam {
TelescopeType GetTelescopeType(const casacore::MeasurementSet& ms) {
  // Read Telescope name and convert to enum
  casacore::ScalarColumn<casacore::String> telescope_name_col(ms.observation(),
                                                              "TELESCOPE_NAME");
  std::string telescope_name = telescope_name_col(0);
  std::for_each(telescope_name.begin(), telescope_name.end(),
                [](char& c) { c = ::toupper(c); });

  if (telescope_name == "AARTFAAC") {
    return kAARTFAAC;
  } else if (telescope_name.compare(0, 4, "ATCA") == 0) {
    return kATCATelescope;
  } else if (telescope_name == "ALMA") {
    return kALMATelescope;
  } else if (telescope_name.compare(0, 4, "EVLA") == 0) {
    return kVLATelescope;
  } else if (telescope_name == "GMRT") {
    return kGMRTTelescope;
  } else if (telescope_name == "LOFAR") {
    return kLofarTelescope;
  } else if (telescope_name == "MEERKAT") {
    return kMeerKATTelescope;
  } else if (telescope_name == "MID") {
    return kSkaMidTelescope;
  } else if (telescope_name == "MWA") {
    return kMWATelescope;
    // check if telescope_name starts with "OSKAR"
  } else if (telescope_name.rfind("OSKAR", 0) == 0) {
    return kOSKARTelescope;
  } else if (telescope_name == "OVRO_MMA" || telescope_name == "OVRO_LWA") {
    return kOvroLwaTelescope;
  } else {
    return kUnknownTelescope;
  }
}

std::unique_ptr<telescope::Telescope> Load(const casacore::MeasurementSet& ms,
                                           const Options& options) {
  std::unique_ptr<telescope::Telescope> telescope;
  const TelescopeType telescope_name = GetTelescopeType(ms);
  switch (telescope_name) {
    case kAARTFAAC:
    case kLofarTelescope:
      telescope = std::make_unique<telescope::LOFAR>(ms, options);
      break;
    case kALMATelescope:
      telescope = std::make_unique<telescope::Alma>(ms, options);
      break;
    case kATCATelescope: {
      auto coefs = std::make_unique<circularsymmetric::ATCACoefficients>();
      telescope =
          std::make_unique<telescope::Dish>(ms, std::move(coefs), options);
    } break;
    case kGMRTTelescope: {
      auto coefs = std::make_unique<circularsymmetric::GMRTCoefficients>();
      telescope =
          std::make_unique<telescope::Dish>(ms, std::move(coefs), options);
    } break;
    case kMeerKATTelescope: {
      auto coefs = std::make_unique<circularsymmetric::MeerKATCoefficients>();
      telescope =
          std::make_unique<telescope::Dish>(ms, std::move(coefs), options);
    } break;
    case kMWATelescope:
      telescope = std::make_unique<telescope::MWA>(ms, options);
      break;
    case kOSKARTelescope:
      telescope = std::make_unique<telescope::OSKAR>(ms, options);
      break;
    case kSkaMidTelescope:
      telescope = std::make_unique<telescope::SkaMid>(ms, options);
      break;
    case kVLATelescope: {
      auto coefs = std::make_unique<circularsymmetric::VLACoefficients>("");
      telescope =
          std::make_unique<telescope::Dish>(ms, std::move(coefs), options);
    } break;
    case kOvroLwaTelescope: {
      telescope = std::make_unique<telescope::Lwa>(ms, options);
    } break;
    default:
      casacore::ScalarColumn<casacore::String> telescope_name_col(
          ms.observation(), "TELESCOPE_NAME");
      std::stringstream message;
      message << "The requested telescope type " << telescope_name_col(0)
              << " is not implemented.";
      throw std::runtime_error(message.str());
  }
  return telescope;
}

std::unique_ptr<telescope::Telescope> Load(const std::string& ms_name,
                                           const Options& options) {
  casacore::MeasurementSet ms(ms_name);
  return Load(ms, options);
}

ElementResponseModel GetElementResponseEnum(
    const std::string& element_response) {
  return ElementResponseModelFromString(element_response);
}
}  // namespace everybeam
