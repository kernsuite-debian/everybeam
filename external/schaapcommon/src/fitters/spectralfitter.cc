// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "spectralfitter.h"

#include <limits>

#include "polynomialfitter.h"
#include "nlplfitter.h"

using aocommon::Image;

namespace schaapcommon {
namespace fitters {

namespace {
constexpr double kDefaultReferenceFrequency = 150.0e6;
}

void SpectralFitter::SetFrequencies(const double* frequencies,
                                    const NumT* weights, size_t n) {
  frequencies_.assign(frequencies, frequencies + n);
  weights_.assign(weights, weights + n);
  NumT weight_sum = 0.0;
  reference_frequency_ = 0.0;
  for (size_t i = 0; i != n; ++i) {
    reference_frequency_ += frequencies_[i] * weights_[i];
    weight_sum += weights_[i];
  }
  if (weight_sum > 0.0) {
    reference_frequency_ /= weight_sum;
  } else {
    reference_frequency_ = kDefaultReferenceFrequency;
  }
}

void SpectralFitter::Fit(std::vector<NumT>& terms, const NumT* values, size_t x,
                         size_t y) const {
  if (IsForced()) {
    ForcedFit(terms, values, x, y);
  } else {
    switch (mode_) {
      default:
      case SpectralFittingMode::NoFitting:
        break;

      case SpectralFittingMode::Polynomial: {
        PolynomialFitter fitter;
        const double reference = ReferenceFrequency();
        for (size_t i = 0; i != frequencies_.size(); ++i) {
          if (weights_[i] > 0.0) {
            fitter.AddDataPoint(frequencies_[i] / reference - 1.0, values[i],
                                weights_[i]);
          }
        }

        fitter.Fit(terms, n_terms_);
      } break;

      case SpectralFittingMode::LogPolynomial: {
        NonLinearPowerLawFitter fitter;
        const double reference = ReferenceFrequency();
        for (size_t i = 0; i != frequencies_.size(); ++i) {
          if (weights_[i] > 0.0) {
            fitter.AddDataPoint(frequencies_[i] / reference, values[i]);
          }
        }

        fitter.Fit(terms, n_terms_);
      } break;
    }
  }
}

void SpectralFitter::ForcedFit(std::vector<NumT>& terms,
                               const SpectralFitter::NumT* values, size_t x,
                               size_t y) const {
  terms.resize(n_terms_);
  // We need to find alpha such that
  // y[i] = A f(x[i], terms), with f the shape.
  // The least-squares fit is:
  // A = sum (y[i] w[i] f[i]) / sum (w[i] f[i]^2)
  // However, it turns out that finding the true least-squares solution for A
  // leads to unstable cleaning. This is because a LS constrained flux might
  // integrate to zero. If it does, the peak finding that uses integrated
  // values will again find the same peak (over and over...). Therefore,
  // we now use the linear average to estimate the flux:
  // A = sum (y[i] w[i]) / sum (w[i] f[i])
  // This is what is calculated below.
  terms[0] = 1.0;
  for (size_t term = 1; term != n_terms_; ++term) {
    const Image& term_image = forced_terms_[term - 1];
    terms[term] = term_image[x + y * term_image.Width()];
  }
  float numerator = 0.0;
  float divisor = 0.0;
  for (size_t i = 0; i != frequencies_.size(); ++i) {
    const float w = weights_[i];
    const float f = NonLinearPowerLawFitter::Evaluate(frequencies_[i], terms,
                                                      ReferenceFrequency());
    if (w > 0.0) {
      numerator += w * values[i];
      divisor += w * f;
    }
  }
  terms[0] = numerator / divisor;
}

void SpectralFitter::Evaluate(NumT* values,
                              const std::vector<NumT>& terms) const {
  switch (mode_) {
    default:
    case SpectralFittingMode::NoFitting:
      break;

    case SpectralFittingMode::Polynomial: {
      const double reference = ReferenceFrequency();
      for (size_t i = 0; i != frequencies_.size(); ++i) {
        values[i] = PolynomialFitter::Evaluate(
            frequencies_[i] / reference - 1.0, terms);
      }
    } break;

    case SpectralFittingMode::LogPolynomial: {
      const double reference = ReferenceFrequency();
      for (size_t i = 0; i != frequencies_.size(); ++i) {
        values[i] = NonLinearPowerLawFitter::Evaluate(frequencies_[i], terms,
                                                      reference);
      }
    } break;
  }
}

SpectralFitter::NumT SpectralFitter::Evaluate(const std::vector<NumT>& terms,
                                              double frequency) const {
  switch (mode_) {
    default:
    case SpectralFittingMode::NoFitting:
      throw std::runtime_error(
          "Something is inconsistent: can't evaluate terms at frequency "
          "without fitting");

    case SpectralFittingMode::Polynomial:
      return PolynomialFitter::Evaluate(frequency / ReferenceFrequency() - 1.0,
                                        terms);

    case SpectralFittingMode::LogPolynomial:
      return NonLinearPowerLawFitter::Evaluate(frequency, terms,
                                               ReferenceFrequency());
  }
}

}  // namespace fitters
}  // namespace schaapcommon