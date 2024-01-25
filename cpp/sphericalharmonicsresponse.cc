// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "sphericalharmonicsresponse.h"

#include <array>
#include <filesystem>

#include <H5Cpp.h>

#include "common/mathutils.h"
#include "common/sphericalharmonics.h"
#include "sphericalharmonicsresponsefixeddirection.h"

namespace everybeam {

namespace {
H5::H5File OpenH5File(const std::string& coefficients_file) {
  if (!std::filesystem::exists(coefficients_file)) {
    throw std::runtime_error("Coefficients file " + coefficients_file +
                             " does not exist.");
  }

  H5::H5File h5file;
  try {
    h5file.openFile(coefficients_file.c_str(), H5F_ACC_RDONLY);
  } catch (const H5::FileIException& e) {
    throw std::runtime_error("Could not open coefficients file: " +
                             coefficients_file);
  }

  return h5file;
}

H5::DataSet OpenDataSet(H5::H5File& h5file, const std::string& name) {
  try {
    return h5file.openDataSet(name);
  } catch (H5::Exception& exception) {
    throw std::runtime_error("Error reading dataset " + name + " from " +
                             h5file.getFileName() + ": " +
                             exception.getDetailMsg());
  }
}

template <typename... Types>
void ReadDataSet(H5::DataSet& dataset, Types... arguments) {
  try {
    dataset.read(arguments...);
  } catch (H5::Exception& exception) {
    throw std::runtime_error("Error reading HDF5 dataset : " +
                             exception.getDetailMsg());
  }
}

xt::xtensor<std::complex<double>, 4> ReadCoefficients(
    H5::H5File& h5file, std::optional<std::size_t> element_index) {
  // Construct an H5 type that represents a complex double.
  H5::CompType h5_complex_double_type(sizeof(std::complex<double>));
  h5_complex_double_type.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
  h5_complex_double_type.insertMember("i", sizeof(double),
                                      H5::PredType::NATIVE_DOUBLE);

  H5::DataSet dataset = OpenDataSet(h5file, "coefficients");
  H5::DataSpace dataspace = dataset.getSpace();

  const int n_dimensions = dataspace.getSimpleExtentNdims();
  if (n_dimensions != 4) {
    throw std::runtime_error(
        "Coefficient file " + h5file.getFileName() +
        " has invalid number of dimensions: " + std::to_string(n_dimensions));
  }

  std::array<hsize_t, 4> h5shape;
  dataspace.getSimpleExtentDims(h5shape.data(), nullptr);

  xt::xtensor<std::complex<double>, 4> coefficients;
  if (element_index) {
    // Define the part of the coefficients to read.
    if (*element_index >= h5shape[2]) {
      throw std::runtime_error(
          "Element index " + std::to_string(*element_index) +
          " not found in coefficient file " + h5file.getFileName());
    }
    h5shape[2] = 1;
    const std::array<hsize_t, 4> h5offsets{0, 0, *element_index, 0};
    dataspace.selectHyperslab(H5S_SELECT_SET, h5shape.data(), h5offsets.data());

    H5::DataSpace memspace{4, h5shape.data()};

    const std::array<std::size_t, 4> shape{h5shape[0], h5shape[1], 1,
                                           h5shape[3]};
    coefficients.resize(shape);
    ReadDataSet(dataset, coefficients.data(), h5_complex_double_type, memspace,
                dataspace);
  } else {
    const std::array<std::size_t, 4> shape{h5shape[0], h5shape[1], h5shape[2],
                                           h5shape[3]};
    coefficients.resize(shape);
    ReadDataSet(dataset, coefficients.data(), h5_complex_double_type);
  }

  return coefficients;
}

xt::xtensor<double, 1> ReadFrequencies(H5::H5File& h5file) {
  H5::DataSet dataset = OpenDataSet(h5file, "frequencies");
  H5::DataSpace dataspace = dataset.getSpace();

  const int n_dimensions = dataspace.getSimpleExtentNdims();
  if (n_dimensions != 1) {
    throw std::runtime_error(
        "Coefficient file " + h5file.getFileName() +
        " has invalid number of dimensions for frequencies. (" +
        std::to_string(n_dimensions) + " should be 1).");
  }

  const std::size_t n_elements = dataspace.getSimpleExtentNpoints();
  xt::xtensor<double, 1> frequencies{std::array<std::size_t, 1>{n_elements}};
  ReadDataSet(dataset, frequencies.data(), H5::PredType::NATIVE_DOUBLE);

  // Check if 'frequencies' is sorted and does not contain duplicates.
  for (std::size_t i = 1; i < frequencies.size(); ++i) {
    if (frequencies(i - 1) >= frequencies(i)) {
      throw std::runtime_error("Frequencies in coefficient file " +
                               h5file.getFileName() + " are not sorted.");
    }
  }

  return frequencies;
}

xt::xtensor<int, 2> ReadNms(H5::H5File& h5file) {
  H5::DataSet dataset = OpenDataSet(h5file, "nms");
  H5::DataSpace dataspace = dataset.getSpace();

  const int n_dimensions = dataspace.getSimpleExtentNdims();
  if (n_dimensions != 2) {
    throw std::runtime_error("Coefficient file " + h5file.getFileName() +
                             " has invalid number of dimensions " +
                             std::to_string(n_dimensions) + " for NMS.");
  }

  std::array<hsize_t, 2> h5shape;
  dataspace.getSimpleExtentDims(h5shape.data(), nullptr);
  if (h5shape[1] != 3) {
    throw std::runtime_error("Coefficient file " + h5file.getFileName() +
                             " has invalid size " + std::to_string(h5shape[1]) +
                             " for second NMS dimension.");
  }

  xt::xtensor<int, 2> nms{std::array<std::size_t, 2>{h5shape[0], h5shape[1]}};
  ReadDataSet(dataset, nms.data(), H5::PredType::NATIVE_INT);
  return nms;
}

}  // namespace

SphericalHarmonicsResponse::SphericalHarmonicsResponse(
    const std::string& coefficients_file,
    std::optional<std::size_t> element_index)
    : element_index_(element_index) {
  H5::H5File h5file = OpenH5File(coefficients_file);
  coefficients_ = ReadCoefficients(h5file, element_index);
  frequencies_ = ReadFrequencies(h5file);
  nms_ = ReadNms(h5file);
  h5file.close();

  if ((coefficients_.shape(0) != 2) ||
      (coefficients_.shape(1) != frequencies_.shape(0)) ||
      (coefficients_.shape(3) != nms_.shape(0))) {
    throw std::runtime_error("Inconsistent shape(s) in coefficient file " +
                             coefficients_file);
  }
}

aocommon::MC2x2 SphericalHarmonicsResponse::Response(double frequency,
                                                     double theta,
                                                     double phi) const {
  if (!element_index_) {
    throw std::runtime_error(
        "SphericalHarmonicsResponse needs an element id, since it loaded "
        "coefficients for all elements.");
  }
  return ComputeResponse(0, frequency, theta, phi);
}

aocommon::MC2x2 SphericalHarmonicsResponse::Response(int element_id,
                                                     double frequency,
                                                     double theta,
                                                     double phi) const {
  const std::size_t element_index = element_id;
  if (element_index_) {
    if (*element_index_ != element_index) {
      throw std::runtime_error(
          "Requested element " + std::to_string(element_id) +
          " does not match loaded element " + std::to_string(*element_index_));
    }
  } else {
    if (element_index >= coefficients_.shape(2)) {
      throw std::runtime_error("Element id " + std::to_string(element_id) +
                               " is out of range");
    }
  }
  return ComputeResponse(element_index, frequency, theta, phi);
}

aocommon::MC2x2 SphericalHarmonicsResponse::ComputeResponse(
    std::size_t element_index, double frequency, double theta,
    double phi) const {
  const std::size_t frequency_index = FindFrequencyIndex(frequency);
  aocommon::MC2x2 response = aocommon::MC2x2::Zero();

  for (std::size_t i = 0; i < nms_.shape(0); ++i) {
    // TODO: Vectorize using MatrixComplexDouble2x2 when aocommon has it.
    const std::complex<double> c0 =
        coefficients_(0, frequency_index, element_index, i);
    const std::complex<double> c1 =
        coefficients_(1, frequency_index, element_index, i);
    std::complex<double> q2;
    std::complex<double> q3;
    std::tie(q2, q3) =
        common::F4far_new(nms_(i, 2), nms_(i, 1), nms_(i, 0), theta, phi);
    response[0] += q2 * c0;  // xx
    response[1] += q3 * c0;  // xy
    response[2] += q2 * c1;  // yx
    response[3] += q3 * c1;  // yy
  }

  return response;
}

std::size_t SphericalHarmonicsResponse::FindFrequencyIndex(
    double frequency) const {
  // Find the closest frequency.
  // TODO: Exploit that the frequency list is sorted.
  auto is_closer = [frequency](double x, double y) {
    return std::abs(x - frequency) < std::abs(y - frequency);
  };
  auto result =
      std::min_element(frequencies_.begin(), frequencies_.end(), is_closer);
  return std::distance(frequencies_.begin(), result);
}

std::shared_ptr<ElementResponse> SphericalHarmonicsResponse::FixateDirection(
    const vector3r_t& direction) const {
  const vector2r_t thetaphi = cart2thetaphi(direction);

  return std::make_shared<SphericalHarmonicsResponseFixedDirection>(
      std::static_pointer_cast<const SphericalHarmonicsResponse>(
          shared_from_this()),
      thetaphi[0], thetaphi[1]);
}

}  // namespace everybeam
