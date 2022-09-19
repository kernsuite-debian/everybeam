// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>

#include "config.h"
#include "lobeselementresponse.h"

#include "../common/sphericalharmonics.h"

#include <aocommon/throwruntimeerror.h>

#include <H5Cpp.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/utility/string_view.hpp>
#include <boost/optional.hpp>

#if __cplusplus > 201402L
#include <charconv>
#endif
#include <complex>
#include <map>

// There are two main modi for the AARTFAAC telescope, AARTFAAC-6 and
// AARTFAAC-12. To properly use AARTFAAC in LOBEs mode the coefficients of all
// stations need to be available. At the moment of writing only a partial set
// is available. This means only AARTFAAC-6 is tested.
static const std::array<boost::string_view, 12> kAartfaacStationNames{
    // Available
    "CS002LBA", "CS003LBA", "CS004LBA", "CS005LBA", "CS006LBA", "CS007LBA",
    // Currently unavailable
    "CS001LBA", "CS011LBA", "CS013LBA", "CS017LBA", "CS021LBA", "CS032LBA"};

struct AartfaacStation {
  boost::string_view station;
  int element;
};

#if __cplusplus > 201402L
template <class T>
static T ExtractIntegral(boost::string_view string) {
  int value;
  std::from_chars_result result =
      std::from_chars(string.begin(), string.end(), value);
  if (result.ec != std::errc{} || result.ptr != string.end()) {
    aocommon::ThrowRuntimeError("The value '", string,
                                "' can't be converted to a number");
  }
  return value;
}
#else
/** Modelled after std::from_chars_result. */
struct from_chars_result {
  const char* ptr;
  std::errc ec;
};

/**
 * A minimal implementation of std::from_chars.
 *
 * @note No support for floating point values.
 * @note No support for bases other than 10.
 */
template <class T>
static from_chars_result from_chars(const char* first, const char* last,
                                    T& value) {
  if (first == last || *first < '0' || *first > '9') {
    return {first, std::errc::invalid_argument};
  }
  value = 0;
  do {
    value *= 10;
    value += *first - '0';
    ++first;
  } while (first != last && *first >= '0' && *first <= '9');

  return {last, std::errc{}};
}

template <class T>
static int ExtractIntegral(boost::string_view string) {
  T value;
  from_chars_result result = from_chars(string.begin(), string.end(), value);
  if (result.ec != std::errc{} || result.ptr != string.end()) {
    aocommon::ThrowRuntimeError("The value '", string,
                                "' can't be converted to a number");
  }
  return value;
}
#endif
enum class AartfaacElements { kInner, kOuter };

static boost::optional<AartfaacStation> GetAartfaacStation(
    boost::string_view station_name, AartfaacElements elements) {
  if (!boost::starts_with(station_name, "A12_")) {
    return {};
  }

  station_name.remove_prefix(4);
  const int id = ExtractIntegral<int>(station_name);
  const size_t station_id = id / 48;
  const int element_id =
      id % 48 + (elements == AartfaacElements::kInner ? 0 : 48);

  if (station_id >= kAartfaacStationNames.size()) {
    aocommon::ThrowRuntimeError("Aartfaac station id '", station_id,
                                "' is invalid");
  }
  return AartfaacStation{kAartfaacStationNames[station_id], element_id};
}

namespace everybeam {

namespace {
/**
 * @brief Search for LOBES h5 coefficient file
 * on the suggested path \param search_path. Returns
 * an empty string if the file cannot be found.
 *
 * @param search_path Search path
 * @param station_name Station name, as read from MS
 * @return std::string Path to file or empty string if file cannot be found
 */
boost::filesystem::path FindCoeffFile(const std::string& search_path,
                                      boost::string_view station_name) {
  const std::string station_file = "LOBES_" + std::string{station_name} + ".h5";
  return search_path.empty()
             ? boost::filesystem::path(std::string{EVERYBEAM_DATA_DIR} +
                                       std::string{"/lobes"}) /
                   station_file
             : boost::filesystem::path(search_path) / station_file;
}
}  // namespace

static const H5::CompType kH5Dcomplex = [] {
  const std::string REAL("r");
  const std::string IMAG("i");

  H5::CompType h5_dcomplex(sizeof(std::complex<double>));
  h5_dcomplex.insertMember(REAL, 0, H5::PredType::NATIVE_DOUBLE);
  h5_dcomplex.insertMember(IMAG, sizeof(double), H5::PredType::NATIVE_DOUBLE);
  return h5_dcomplex;
}();

static void ReadAllElements(
    Eigen::Tensor<std::complex<double>, 4, Eigen::RowMajor>& coefficients,
    const H5::DataSet& dataset, const std::vector<unsigned int>& shape) {
  coefficients.resize(shape[0], shape[1], shape[2], shape[3]);
  dataset.read(coefficients.data(), kH5Dcomplex);
}

void ReadOneElement(
    Eigen::Tensor<std::complex<double>, 4, Eigen::RowMajor>& coefficients,
    const H5::DataSet& dataset, const std::vector<unsigned>& shape,
    unsigned index) {
  static constexpr size_t kRank = 4;

  // Define the part of the coefficients to read.
  const std::array<hsize_t, kRank> kOffset = {0, 0, index, 0};
  const std::array<hsize_t, kRank> kCount = {shape[0], shape[1], 1, shape[3]};
  static constexpr std::array<hsize_t, kRank> kStride = {1, 1, 1, 1};
  static constexpr std::array<hsize_t, kRank> kBlock = {1, 1, 1, 1};

  H5::DataSpace memspace{kRank, kCount.data()};
  H5::DataSpace dataspace = dataset.getSpace();
  dataspace.selectHyperslab(H5S_SELECT_SET, kCount.data(), kOffset.data(),
                            kStride.data(), kBlock.data());

  // TODO AST-807 The exact mapping between the data-layout of HD5 and
  // Eigen-Tensors needs to be investigated so the elements can be copied more
  // efficiently.  (Ideally they would be directly read in the proper shape.)
  std::vector<std::complex<double>> buffer(shape[0] * shape[1] * 1 * shape[3]);
  dataset.read(buffer.data(), kH5Dcomplex, memspace, dataspace);
  auto iterator = buffer.begin();
  coefficients.resize(shape[0], shape[1], 1, shape[3]);
  for (size_t i = 0; i < shape[0]; ++i) {
    for (size_t j = 0; j < shape[1]; ++j) {
      for (size_t k = 0; k < shape[3]; ++k) {
        coefficients(i, j, 0, k) = *iterator++;
      }
    }
  }
#if 0
  // TODO AST-807 remove this validation code.
  // This validates the data has been read correctly when compared with the
  // original read function.
  Eigen::Tensor<std::complex<double>, 4, Eigen::RowMajor> expected;
  ReadAllElements(expected, dataset, shape);
  for (size_t i = 0; i < shape[0]; ++i) {
    for (size_t j = 0; j < shape[1]; ++j) {
      for (size_t k = 0; k < shape[3]; ++k) {
        if (coefficients(i, j, 0, k) != expected(i, j, index, k)) {
          asm("int3");
        }
      }
    }
  }
#endif
}

LOBESElementResponse::LOBESElementResponse(const std::string& name,
                                           const Options& options) {
  const boost::optional<AartfaacStation> aartfaac_station =
      GetAartfaacStation(name, AartfaacElements::kInner);

  boost::filesystem::path coeff_file_path = FindCoeffFile(
      options.coeff_path, aartfaac_station ? aartfaac_station->station : name);
  H5::H5File h5file;

  if (!boost::filesystem::exists(coeff_file_path)) {
    throw std::runtime_error("LOBES coeffcients file: " +
                             coeff_file_path.string() + " does not exists");
  }

  try {
    h5file.openFile(coeff_file_path.c_str(), H5F_ACC_RDONLY);
  } catch (const H5::FileIException& e) {
    throw std::runtime_error("Could not open LOBES coeffcients file: " +
                             coeff_file_path.string());
  }

  H5::DataSet dataset = h5file.openDataSet("coefficients");
  H5::DataSpace dataspace = dataset.getSpace();
  int nr_elements = dataspace.getSimpleExtentNpoints();

  // Get the number of dimensions in the dataspace.
  int ndims_coefficients = dataspace.getSimpleExtentNdims();

  // Get the dimension size of each dimension in the dataspace and display them.
  std::vector<hsize_t> dims_coefficients(ndims_coefficients);
  dataspace.getSimpleExtentDims(dims_coefficients.data(), nullptr);
  const std::vector<unsigned int> coefficients_shape(dims_coefficients.begin(),
                                                     dims_coefficients.end());

  if (aartfaac_station) {
    std::stringstream sstr;
    ReadOneElement(coefficients_, dataset, coefficients_shape,
                   aartfaac_station->element);
  } else {
    ReadAllElements(coefficients_, dataset, coefficients_shape);
  }

  // Frequencies
  dataset = h5file.openDataSet("frequencies");
  dataspace = dataset.getSpace();
  nr_elements = dataspace.getSimpleExtentNpoints();

  frequencies_.resize(nr_elements);
  dataset.read(frequencies_.data(), H5::PredType::NATIVE_DOUBLE);

  // nms
  dataset = h5file.openDataSet("nms");
  dataspace = dataset.getSpace();
  nr_elements = dataspace.getSimpleExtentNpoints();

  // Get the number of dimensions in the dataspace.
  int ndims_nms = dataspace.getSimpleExtentNdims();

  // Get the dimension size of each dimension in the dataspace and display them.
  std::vector<hsize_t> dims_nms(ndims_nms);
  dataspace.getSimpleExtentDims(dims_nms.data(), nullptr);

  nms_.resize(dims_nms[0]);
  dataset.read(nms_.data(), H5::PredType::NATIVE_INT);
}

LOBESElementResponse::BaseFunctions LOBESElementResponse::ComputeBaseFunctions(
    double theta, double phi) const {
  LOBESElementResponse::BaseFunctions base_functions(nms_.size(), 2);
  base_functions.setZero();

  for (size_t i = 0; i < nms_.size(); ++i) {
    auto nms = nms_[i];
    std::complex<double> q2, q3;
    std::tie(q2, q3) =
        everybeam::common::F4far_new(nms.s, nms.m, nms.n, theta, phi);
    base_functions(i, 0) = q2;
    base_functions(i, 1) = q3;
  }
  return base_functions;
}

aocommon::MC2x2 LOBESElementResponse::Response(int element_id, double freq,
                                               double theta, double phi) const {
  // Clip directions below the horizon.
  if (theta >= M_PI_2) {
    return aocommon::MC2x2::Zero();
  }

  // When the objects basefunctions_ aren't initialized create our own copy.
  // Note it's not possible to set the object's version since the function is
  // called from multiple threads.
  const BaseFunctions& basefunctions =
      basefunctions_ ? *basefunctions_ : ComputeBaseFunctions(theta, phi);

  const int freq_idx = FindFrequencyIdx(freq);
  std::complex<double> xx = {0}, xy = {0}, yx = {0}, yy = {0};

  const int nr_rows = basefunctions.rows();
  if (nr_rows == 0) {
    throw std::runtime_error(
        "Number of rows in basefunctions_ member is 0. Did you run "
        "SetFieldQuantities?");
  }

  for (int i = 0; i < nr_rows; ++i) {
    const std::complex<double> q2 = basefunctions(i, 0);
    const std::complex<double> q3 = basefunctions(i, 1);
    xx += q2 * coefficients_(0, freq_idx, element_id, i);
    xy += q3 * coefficients_(0, freq_idx, element_id, i);
    yx += q2 * coefficients_(1, freq_idx, element_id, i);
    yy += q3 * coefficients_(1, freq_idx, element_id, i);
  }

  return aocommon::MC2x2(xx, xy, yx, yy);
}

std::shared_ptr<LOBESElementResponse> LOBESElementResponse::GetInstance(
    const std::string& name, const Options& options) {
  static std::map<std::string, std::shared_ptr<LOBESElementResponse>>
      name_response_map;

  auto entry = name_response_map.find(name);
  if (entry == name_response_map.end()) {
    entry = name_response_map.insert(
        entry, {name, std::make_shared<LOBESElementResponse>(name, options)});
  }
  return entry->second;
}

}  // namespace everybeam
