// Options.h: Class for specifying the telescope response options.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_OPTIONS_H_
#define EVERYBEAM_OPTIONS_H_

#include <string>
#include <vector>
#include "elementresponse.h"

namespace everybeam {

/**
 * @brief Class/Struct specifying everybeam Options. Needs further
 * implementation!
 *
 */
struct Options {
  // Path to coefficients file
  std::string coeff_path = ".";

  // LOFAR specific
  bool use_differential_beam = false;
  bool use_channel_frequency = true;
  std::string data_column_name = "DATA";
  ElementResponseModel element_response_model = ElementResponseModel::kHamaker;

  // MWA specific (Lofar probably will follow)
  bool frequency_interpolation = false;
};

struct ATermSettings {
  // Path to coefficients file
  std::string coeff_path = ".";
  // Save aterm fits files?
  bool save_aterms = false;
  // Prefix for the aterm fits files
  std::string save_aterms_prefix = "wsclean";
  // Default for the data column name
  std::string data_column_name = "DATA";

  std::vector<std::string> filenames = std::vector<std::string>();
  size_t max_support = 32;

  // Time interval (s) before a new aterm will be computed
  double aterm_update_interval = 300.;
  size_t padded_image_width = 0, padded_image_height = 0,
         trimmed_image_width = 0, trimmed_image_height = 0;
};
}  // namespace everybeam
#endif  // EVERYBEAM_OPTIONS_H_