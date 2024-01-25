// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oskar.h"

#include <cassert>

#include <aocommon/banddata.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

#include "../common/mathutils.h"
#include "../common/casautils.h"
#include "../msreadutils.h"

using casacore::MeasurementSet;
using everybeam::Station;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::pointresponse::PointResponse;
using everybeam::telescope::OSKAR;

OSKAR::OSKAR(const MeasurementSet& ms, const Options& options)
    : PhasedArray(ms, options) {
  if (GetOptions().element_response_model == ElementResponseModel::kDefault) {
    options_.element_response_model = ElementResponseModel::kOSKARSphericalWave;
  }
  // OSKAR never uses the subband frequency.
  options_.use_channel_frequency = true;

  ReadAllStations(ms, stations_.begin(), GetOptions());

  aocommon::BandData band(ms.spectralWindow());
  casacore::ScalarMeasColumn<casacore::MDirection> delay_dir_col(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));

  CorrectionMode preapplied_correction_mode;
  casacore::MDirection preapplied_beam_dir;
  PhasedArray::CalculatePreappliedBeamOptions(ms, options_.data_column_name,
                                              preapplied_beam_dir,
                                              preapplied_correction_mode);

  size_t channel_count = band.ChannelCount();
  std::vector<double> channel_freqs(channel_count);
  for (size_t idx = 0; idx < channel_count; ++idx) {
    channel_freqs[idx] = band.ChannelFrequency(idx);
  }

  // Populate struct
  ms_properties_ = MSProperties();
  ms_properties_.subband_freq = 0.0;  // Since use_channel_frequency == true
  ms_properties_.delay_dir = delay_dir_col(0);
  // tile_beam_dir has dummy values for OSKAR
  ms_properties_.tile_beam_dir = delay_dir_col(0);
  ms_properties_.preapplied_beam_dir = preapplied_beam_dir;
  ms_properties_.preapplied_correction_mode = preapplied_correction_mode;
  ms_properties_.channel_count = channel_count;
  ms_properties_.channel_freqs = channel_freqs;
}
