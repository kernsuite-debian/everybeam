#ifndef EVERYBEAM_POINTRESPONSE_AIRY_POINT_H_
#define EVERYBEAM_POINTRESPONSE_AIRY_POINT_H_

#include "pointresponse.h"

namespace everybeam {
namespace pointresponse {

/**
 * @brief Class for computing the directional response of telescopes with
 * an Airy Disc response, e.g. ALMA.
 */
class [[gnu::visibility("default")]] AiryPoint final : public PointResponse {
 public:
  AiryPoint(const telescope::Telescope* telescope_ptr, double time)
      : PointResponse(telescope_ptr, time){};

  void Response(BeamMode beam_mode, std::complex<float> * buffer, double ra,
                double dec, double freq, size_t station_idx, size_t field_id)
      override;

  void ResponseAllStations(BeamMode beam_mode, std::complex<float> * buffer,
                           double ra, double dec, double freq, size_t field_id)
      override;
};
}  // namespace pointresponse
}  // namespace everybeam
#endif  // EVERYBEAM_POINTRESPONSE_DISHPOINT_H_
