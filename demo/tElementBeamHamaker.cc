// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "config.h"

#include "tElementBeamCommon.h"

using namespace everybeam;

int main(int argc, char** argv) {
  ElementResponseModel model(kHamaker);
  double frequency = 132e6;  // Mhz
  std::string input_filename(TEST_MEASUREMENTSET);
  std::string output_filename("element-beams-hamaker.fits");
  run(model, frequency, input_filename, output_filename);
}