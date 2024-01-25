// lwa.h: Base class for computing the response for the OVRO-LWA
// telescope.
//
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_LWA_H_
#define EVERYBEAM_TELESCOPE_LWA_H_

#include <aocommon/banddata.h>

#include "phasedarray.h"
#include "../msreadutils.h"

namespace everybeam {

namespace telescope {

// LWA telescope class
class [[gnu::visibility("default")]] Lwa final : public PhasedArray {
 public:
  /**
   * @brief Construct a new Lwa object
   *
   * @param ms MeasurementSet
   * @param options telescope options
   */
  Lwa(const casacore::MeasurementSet& ms, const Options& options)
      : PhasedArray(ms, options) {
    ReadAllStations(ms, stations_.begin(), options_);

    // add MS properties as well
    ms_properties_ = MSProperties();

    // Read frequency information
    aocommon::BandData band(ms.spectralWindow());
    size_t channel_count = band.ChannelCount();
    std::vector<double> channel_freqs(channel_count);
    for (size_t idx = 0; idx < channel_count; ++idx) {
      channel_freqs[idx] = band.ChannelFrequency(idx);
    }

    // Read Field information
    casacore::ScalarMeasColumn<casacore::MDirection> delay_dir_col(
        ms.field(),
        casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));

    casacore::ScalarMeasColumn<casacore::MDirection> reference_dir_col(
        ms.field(),
        casacore::MSField::columnName(casacore::MSFieldEnums::REFERENCE_DIR));

    CorrectionMode preapplied_correction_mode;
    casacore::MDirection preapplied_beam_dir;

    PhasedArray::CalculatePreappliedBeamOptions(ms, options_.data_column_name,
                                                preapplied_beam_dir,
                                                preapplied_correction_mode);

    ms_properties_.subband_freq = band.ReferenceFrequency();
    ms_properties_.delay_dir = delay_dir_col(0);
    ms_properties_.reference_dir = reference_dir_col(0);
    ms_properties_.preapplied_beam_dir = preapplied_beam_dir;
    ms_properties_.preapplied_correction_mode = preapplied_correction_mode;
    ms_properties_.channel_count = channel_count;
    ms_properties_.channel_freqs = channel_freqs;
  };
};

}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_LWA_H_
