// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "airypoint.h"
#include "../telescope/alma.h"
#include "../circularsymmetric/voltagepattern.h"
#include "../circularsymmetric/vlacoefficients.h"

#include <aocommon/uvector.h>

namespace everybeam::pointresponse {

void AiryPoint::Response(BeamMode /* beam_mode */, std::complex<float>* buffer,
                         double ra, double dec, double frequency,
                         size_t station_idx, size_t field_id) {
  const telescope::Alma& alma_telescope =
      static_cast<const telescope::Alma&>(GetTelescope());

  double pdir_ra;
  double pdir_dec;
  std::tie(pdir_ra, pdir_dec) = alma_telescope.GetFieldPointing()[field_id];
  const telescope::AiryParameters parameters =
      alma_telescope.GetAiryParameters(station_idx);
  circularsymmetric::VoltagePattern vp({frequency},
                                       parameters.maximum_radius_arc_min);
  vp.EvaluateAiryDisk(parameters.dish_diameter_in_m,
                      parameters.blocked_diameter_in_m);
  vp.Render(buffer, ra, dec, pdir_ra, pdir_dec, frequency);
}

void AiryPoint::ResponseAllStations(BeamMode beam_mode,
                                    std::complex<float>* buffer, double ra,
                                    double dec, double freq, size_t field_id) {
  const telescope::Alma& alma_telescope =
      static_cast<const telescope::Alma&>(GetTelescope());
  if (alma_telescope.IsHomogeneous())
    HomogeneousAllStationsResponse(beam_mode, buffer, ra, dec, freq, field_id);
  else
    ResponseAllStations(beam_mode, buffer, ra, dec, freq, field_id);
}

}  // namespace everybeam::pointresponse
