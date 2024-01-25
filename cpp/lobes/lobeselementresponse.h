// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_LOBES_ELEMENTRESPONSE_H_
#define EVERYBEAM_LOBES_ELEMENTRESPONSE_H_

#include "../options.h"
#include "../sphericalharmonicsresponse.h"

namespace everybeam {

//! Implementation of the Lobes response model
class LOBESElementResponse : public SphericalHarmonicsResponse {
 public:
  /**
   * @brief Construct a new LOBESElementResponse object
   *
   * @param name (LOFAR) station name, i.e. CS302LBA
   * @param options if options.coeff_path is non-empty it is used to find
   * coefficient files
   */
  LOBESElementResponse(const std::string& name, const Options& options)
      : SphericalHarmonicsResponse(SphericalHarmonicsArguments(name, options)) {
  }

  ElementResponseModel GetModel() const final override {
    return ElementResponseModel::kLOBES;
  }

  /**
   * @brief Create LOBESElementResponse
   *
   * @param name Station name, e.g. CS302LBA
   */
  [[gnu::visibility(
      "default")]] static std::shared_ptr<const LOBESElementResponse>
  GetInstance(const std::string& name, const Options& options);

 private:
  static std::tuple<std::string, std::optional<std::size_t>>
  SphericalHarmonicsArguments(const std::string& name, const Options& options);
};

}  // namespace everybeam

#endif
