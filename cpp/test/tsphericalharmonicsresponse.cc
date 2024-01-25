// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// This file tests both SphericalHarmonicsResponse and
// SphericalHarmonicsResponseFixedDirection since they are closely related.
#include "../sphericalharmonicsresponse.h"
#include "../sphericalharmonicsresponsefixeddirection.h"

#include <filesystem>
#include <fstream>
#include <type_traits>

#include <boost/test/unit_test.hpp>
#include <xtensor/xview.hpp>

#include <H5Cpp.h>

#include "common/mathutils.h"

using everybeam::ElementResponse;
using everybeam::ElementResponseModel;
using everybeam::SphericalHarmonicsResponse;

namespace {

const std::string kCoefficientsName{"coefficients"};
const std::string kFrequenciesName{"frequencies"};
const std::string kNmsName{"nms"};

/// Derive a class from SphericalHarmonicsResponse, since it's abstract.
class TestResponse : public SphericalHarmonicsResponse {
 public:
  explicit TestResponse(const std::string& coefficients_file,
                        std::optional<std::size_t> element_index)
      : SphericalHarmonicsResponse(coefficients_file, element_index) {}

  ElementResponseModel GetModel() const override {
    if (model_ == ElementResponseModel::kDefault) {
      BOOST_FAIL("Unexpected GetModel() call");
    }
    return model_;
  }

  void SetModel(ElementResponseModel model) { model_ = model; };

 private:
  ElementResponseModel model_ = ElementResponseModel::kDefault;
};

/// Fixture with static functions for testing if SphericalHarmonicsResponse
/// reads HDF5 files with coefficients etc. correctly.
class H5Fixture {
 public:
  H5Fixture() {}

  ~H5Fixture() {
    for (const std::string& filename : created_files_) {
      std::filesystem::remove(filename);
    }
  }

  template <typename ContentType>
  static void AddDataSet(H5::H5File& h5file, const std::string& name,
                         const ContentType& content) {
    H5::DataType datatype;
    if (std::is_same_v<typename ContentType::value_type,
                       std::complex<double>>) {
      H5::CompType complex_double_type(sizeof(std::complex<double>));
      complex_double_type.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
      complex_double_type.insertMember("i", sizeof(double),
                                       H5::PredType::NATIVE_DOUBLE);
      datatype = complex_double_type;
    } else if (std::is_same_v<typename ContentType::value_type, double>) {
      datatype = H5::PredType::NATIVE_DOUBLE;
    } else if (std::is_same_v<typename ContentType::value_type, int>) {
      datatype = H5::PredType::NATIVE_INT;
    } else {
      BOOST_FAIL("Unsupported data type in test");
    }

    std::array<hsize_t, ContentType::rank> h5shape;
    for (std::size_t dimension = 0; dimension < ContentType::rank;
         ++dimension) {
      h5shape[dimension] = content.shape(dimension);
    }
    H5::DataSpace dataspace(ContentType::rank, h5shape.data());

    H5::DataSet dataset =
        h5file.createDataSet(name.c_str(), datatype, dataspace);
    dataset.write(content.data(), datatype);
  }

  template <typename CoefficientsType, typename FrequenciesType,
            typename NmsType>
  void CreateH5File(const std::string& filename,
                    const CoefficientsType& coefficients,
                    const FrequenciesType& frequencies, const NmsType& nms) {
    H5::H5File h5file(filename, H5F_ACC_TRUNC);
    AddDataSet(h5file, kCoefficientsName, coefficients);
    AddDataSet(h5file, kFrequenciesName, frequencies);
    AddDataSet(h5file, kNmsName, nms);
    h5file.close();
    created_files_.push_back(filename);
  }

  template <typename CoefficientsType, typename FrequenciesType,
            typename NmsType>
  void CheckValid(const CoefficientsType& coefficients,
                  const FrequenciesType& frequencies, const NmsType& nms,
                  std::optional<std::size_t> element_index = std::nullopt) {
    const std::string kFilename = "test-elemresp-spherical-wave-valid.h5";
    CreateH5File(kFilename, coefficients, frequencies, nms);

    const TestResponse element_response{kFilename, element_index};
    BOOST_TEST(element_response.HasFixedElementIndex() ==
               static_cast<bool>(element_index));
    if (!element_index) {
      BOOST_TEST(element_response.GetCoefficients() == coefficients);
    } else {
      // Using xt::range ensures that the number of dimensions remains equal.
      const CoefficientsType expected_coefficients =
          xt::view(coefficients, xt::all(), xt::all(),
                   xt::range(*element_index, *element_index + 1), xt::all());
      BOOST_TEST(element_response.GetCoefficients() == expected_coefficients);
      BOOST_TEST(element_response.GetElementIndex() == *element_index);
    }
    BOOST_TEST(element_response.GetFrequencies() == frequencies);
    BOOST_TEST(element_response.GetNms() == nms);
  }

  template <typename CoefficientsType, typename FrequenciesType,
            typename NmsType>
  void CheckInvalid(const CoefficientsType& coefficients,
                    const FrequenciesType& frequencies, const NmsType& nms,
                    std::optional<std::size_t> element_index = std::nullopt) {
    const std::string kFilename = "test-elemresp-spherical-wave-invalid.h5";
    H5::H5File h5file(kFilename, H5F_ACC_TRUNC);
    if (coefficients.size() > 0)
      AddDataSet(h5file, kCoefficientsName, coefficients);
    if (frequencies.size() > 0)
      AddDataSet(h5file, kFrequenciesName, frequencies);
    if (nms.size() > 0) AddDataSet(h5file, kNmsName, nms);
    h5file.close();
    created_files_.push_back(kFilename);

    BOOST_CHECK_THROW(TestResponse(kFilename, element_index),
                      std::runtime_error);
  }

  static void CheckEqualResponse(const aocommon::MC2x2& left,
                                 const aocommon::MC2x2& right) {
    for (std::size_t i = 0; i < 4; ++i) {
      BOOST_CHECK_CLOSE(left[i].real(), right[i].real(), 1.0e-9);
      BOOST_CHECK_CLOSE(left[i].imag(), right[i].imag(), 1.0e-9);
    }
  }

 private:
  std::vector<std::string> created_files_;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(element_response_spherical_wave)

BOOST_FIXTURE_TEST_CASE(read_coefficients, H5Fixture) {
  // Coefficients should have shape (2, A, B, C)
  const xt::xtensor<std::complex<double>, 4> kValidCoefficients{
      {{{std::complex<double>{5.0, 6.0}}}},
      {{{std::complex<double>{7.0, 8.0}}}},
  };
  const xt::xtensor<std::complex<double>, 4> kEmptyCoefficients;
  const xt::xtensor<std::complex<double>, 5> k5DimensionalCoefficients{
      {{{{std::complex<double>{5.0, 6.0}}}}},
      {{{{std::complex<double>{7.0, 8.0}}}}},
  };
  const xt::xtensor<std::complex<double>, 4> kFirstDimensionInvalidCoefficients{
      {{{std::complex<double>()}}}};
  const xt::xtensor<int, 4> kInvalidTypeCoefficients{{{{5}}}, {{{6}}}};

  // Frequencies should have size A.
  const xt::xtensor<double, 1> kFrequencies{42.0e6};
  // NMS should have size (C, 3).
  const xt::xtensor<int, 2> kNms{{16, 17, 18}};

  CheckValid(kValidCoefficients, kFrequencies, kNms);
  CheckInvalid(kEmptyCoefficients, kFrequencies, kNms);
  CheckInvalid(k5DimensionalCoefficients, kFrequencies, kNms);
  CheckInvalid(kFirstDimensionInvalidCoefficients, kFrequencies, kNms);
  CheckInvalid(kInvalidTypeCoefficients, kFrequencies, kNms);
}

BOOST_FIXTURE_TEST_CASE(read_single_element, H5Fixture) {
  const xt::xtensor<std::complex<double>, 4> kCoefficients{
      {{{std::complex<double>{11.0, 12.0}},
        {std::complex<double>{13.0, 14.0}},
        {std::complex<double>{15.0, 16.0}}}},
      {{{std::complex<double>{21.0, 22.0}},
        {std::complex<double>{23.0, 24.0}},
        {std::complex<double>{25.0, 26.0}}}}};
  const xt::xtensor<double, 1> kFrequencies{42.0e6};
  const xt::xtensor<int, 2> kNms{{16, 17, 18}};

  const std::array<std::size_t, 4> kCoefficientsShape{2, 1, 3, 1};
  BOOST_TEST(kCoefficients.shape() == kCoefficientsShape);  // sanity check

  for (std::size_t element_index = 0; element_index < 3; ++element_index) {
    CheckValid(kCoefficients, kFrequencies, kNms, element_index);
  }
  CheckInvalid(kCoefficients, kFrequencies, kNms, 3);
}

BOOST_FIXTURE_TEST_CASE(read_frequencies, H5Fixture) {
  const xt::xtensor<std::complex<double>, 4> kCoefficients{
      {{{std::complex<double>{2.0, 3.0}}},
       {{std::complex<double>{4.0, 5.0}}},
       {{std::complex<double>{5.0, 6.0}}}},
      {{{std::complex<double>{7.0, 8.0}}},
       {{std::complex<double>{9.0, 10.0}}},
       {{std::complex<double>{11.0, 12.0}}}},
  };
  const xt::xtensor<int, 2> kNms{{16, 17, 18}};

  const xt::xtensor<double, 1> kValidFrequencies{42.0e6, 43.0e6, 44.0e6};
  const xt::xtensor<double, 1> kEmptyFrequencies;
  const xt::xtensor<double, 2> kTwoDimensionalFrequencies{
      {42.0e6, 43.0e6, 44.0e6}};
  const xt::xtensor<double, 1> kInvalidSizeFrequencies{42.0e6, 43.0e6};
  const xt::xtensor<double, 1> kUnsortedFrequencies{43.0e6, 42.0e6, 41.0e6};
  const xt::xtensor<double, 1> kPartiallySortedFrequencies{40.0e6, 42.0e6,
                                                           41.0e6};
  const xt::xtensor<double, 1> kNonUniqueFrequencies{41.0e6, 42.0e6, 42.0e6};
  const xt::xtensor<std::complex<double>, 1> kInvalidTypeFrequencies{
      std::complex<double>{42e6}};

  CheckValid(kCoefficients, kValidFrequencies, kNms);
  CheckInvalid(kCoefficients, kEmptyFrequencies, kNms);
  CheckInvalid(kCoefficients, kTwoDimensionalFrequencies, kNms);
  CheckInvalid(kCoefficients, kInvalidSizeFrequencies, kNms);
  CheckInvalid(kCoefficients, kUnsortedFrequencies, kNms);
  CheckInvalid(kCoefficients, kPartiallySortedFrequencies, kNms);
  CheckInvalid(kCoefficients, kNonUniqueFrequencies, kNms);
  CheckInvalid(kCoefficients, kInvalidTypeFrequencies, kNms);
}

BOOST_FIXTURE_TEST_CASE(read_nms, H5Fixture) {
  const xt::xtensor<std::complex<double>, 4> kCoefficients{
      {{{std::complex<double>{5.0, 6.0}}}},
      {{{std::complex<double>{7.0, 8.0}}}},
  };
  const xt::xtensor<double, 1> kFrequencies{42.0e6};

  const xt::xtensor<int, 2> kValidNms{{16, 17, 18}};
  const xt::xtensor<int, 2> kEmptyNms;
  const xt::xtensor<int, 3> kThreeDimensionalNms{{{1, 2, 3}, {4, 5, 6}}};
  const xt::xtensor<int, 2> kInvalidFirstDimensionNms{{1, 2, 3}, {4, 5, 6}};
  const xt::xtensor<int, 2> kInvalidSecondDimensionNms{{16, 17, 18, 19}};
  const xt::xtensor<std::complex<double>, 2> kInvalidTypeNms{{16, 17, 18}};

  CheckValid(kCoefficients, kFrequencies, kValidNms);
  CheckInvalid(kCoefficients, kFrequencies, kEmptyNms);
  CheckInvalid(kCoefficients, kFrequencies, kThreeDimensionalNms);
  CheckInvalid(kCoefficients, kFrequencies, kInvalidFirstDimensionNms);
  CheckInvalid(kCoefficients, kFrequencies, kInvalidSecondDimensionNms);
  CheckInvalid(kCoefficients, kFrequencies, kInvalidTypeNms);
}

BOOST_AUTO_TEST_CASE(missing_file) {
  BOOST_CHECK_THROW(TestResponse("does_not_exist.h5", std::nullopt),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(invalid_file) {
  const std::string kFilename{"invalid_contents.h5"};
  std::filesystem::path path{kFilename};
  std::ofstream(path) << "This file is not a valid HDF5 file!";
  BOOST_CHECK_THROW(TestResponse(kFilename, std::nullopt), std::runtime_error);
  std::filesystem::remove(path);
}

BOOST_FIXTURE_TEST_CASE(find_frequency_index_multiple, H5Fixture) {
  const std::array<std::size_t, 4> kCoefficientsShape{2, 3, 1, 1};
  const xt::xtensor<std::complex<double>, 4> kCoefficients{kCoefficientsShape};
  const xt::xtensor<int, 2> kNms{{16, 17, 18}};
  const xt::xtensor<double, 1> kFrequencies{42.0e6, 43.0e6, 44.0e6};

  const std::string kFilename = "test_find_freq_index_multiple.h5";
  CreateH5File(kFilename, kCoefficients, kFrequencies, kNms);

  const TestResponse element_response(kFilename, std::nullopt);
  BOOST_TEST(element_response.FindFrequencyIndex(0.0) == 0);
  BOOST_TEST(element_response.FindFrequencyIndex(42.0e6) == 0);
  BOOST_TEST(element_response.FindFrequencyIndex(42.4e6) == 0);
  BOOST_TEST(element_response.FindFrequencyIndex(42.6e6) == 1);
  BOOST_TEST(element_response.FindFrequencyIndex(42.6e6) == 1);
  BOOST_TEST(element_response.FindFrequencyIndex(43.0e6) == 1);
  BOOST_TEST(element_response.FindFrequencyIndex(44.0e6) == 2);
  BOOST_TEST(element_response.FindFrequencyIndex(100.0e6) == 2);
}

BOOST_FIXTURE_TEST_CASE(find_frequency_index_single, H5Fixture) {
  const std::array<std::size_t, 4> kCoefficientsShape{2, 1, 1, 1};
  const xt::xtensor<std::complex<double>, 4> kCoefficients{kCoefficientsShape};
  const xt::xtensor<int, 2> kNms{{16, 17, 18}};
  const xt::xtensor<double, 1> kFrequencies{42.0e6};

  const std::string kFilename = "test_find_freq_index_single.h5";
  CreateH5File(kFilename, kCoefficients, kFrequencies, kNms);

  const TestResponse element_response(kFilename, std::nullopt);
  BOOST_TEST(element_response.FindFrequencyIndex(0.0) == 0);
  BOOST_TEST(element_response.FindFrequencyIndex(42.0e6) == 0);
  BOOST_TEST(element_response.FindFrequencyIndex(100.0e6) == 0);
}

BOOST_FIXTURE_TEST_CASE(fixate_direction, H5Fixture) {
  // Tests the FixateDirection calls of both SphericalHarmonicsResponse and
  // SphericalHarmonicsResponseFixedDirection.

  const xt::xtensor<std::complex<double>, 4> kCoefficients{
      {{{std::complex<double>{1.0, -1.0}}}},
      {{{std::complex<double>{2.0, -2.0}}}},
  };
  const double kFrequency{42.0e6};
  const xt::xtensor<double, 1> kFrequencies{kFrequency};
  const xt::xtensor<int, 2> kNms{{1, -1, 2}};
  const std::size_t kElementIndex = 0;

  const std::string kFilename = "test_fix_direction.h5";
  CreateH5File(kFilename, kCoefficients, kFrequencies, kNms);
  const auto element_response =
      std::make_shared<TestResponse>(kFilename, kElementIndex);

  const std::array<double, 3> kDirection1{20.0, -10.0, 5.0};
  const std::array<double, 3> kDirection2{3.0, 42.0, -1.0};
  const std::array<double, 2> kThetaPhi1 =
      everybeam::cart2thetaphi(kDirection1);
  const std::array<double, 2> kThetaPhi2 =
      everybeam::cart2thetaphi(kDirection2);
  const double kTheta1 = kThetaPhi1[0];
  const double kPhi1 = kThetaPhi1[1];
  const double kTheta2 = kThetaPhi2[0];
  const double kPhi2 = kThetaPhi2[1];

  std::shared_ptr<ElementResponse> fixed_direction_1 =
      element_response->FixateDirection(kDirection1);
  BOOST_TEST(dynamic_cast<everybeam::SphericalHarmonicsResponseFixedDirection*>(
      fixed_direction_1.get()));

  const aocommon::MC2x2 expected_response_1 =
      element_response->Response(kFrequency, kTheta1, kPhi1);
  CheckEqualResponse(expected_response_1,
                     fixed_direction_1->Response(kFrequency, kTheta1, kPhi1));
  CheckEqualResponse(expected_response_1,
                     fixed_direction_1->Response(kFrequency, kTheta2, kPhi2));

  // Test calling FixateDirection on fixed_direction_1.
  std::shared_ptr<ElementResponse> fixed_direction_2 =
      fixed_direction_1->FixateDirection(kDirection2);
  BOOST_TEST(dynamic_cast<everybeam::SphericalHarmonicsResponseFixedDirection*>(
      fixed_direction_2.get()));

  const aocommon::MC2x2 expected_response_2 =
      element_response->Response(kFrequency, kTheta2, kPhi2);
  // Verify that the response is completely different for direction2.
  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_TEST(std::abs(expected_response_1[i].real() -
                        expected_response_2[i].real()) > 0.1);
    BOOST_TEST(std::abs(expected_response_1[i].imag() -
                        expected_response_2[i].imag()) > 0.1);
  }

  CheckEqualResponse(expected_response_2,
                     fixed_direction_2->Response(kFrequency, kTheta1, kPhi1));
  CheckEqualResponse(expected_response_2,
                     fixed_direction_2->Response(kFrequency, kTheta2, kPhi2));
}

BOOST_FIXTURE_TEST_CASE(get_model, H5Fixture) {
  // Test that SphericalHarmonicsResponseFixedDirection::GetModel calls
  // GetModel of the class that generated it.

  const xt::xtensor<std::complex<double>, 4> kCoefficients{
      {{{std::complex<double>{1.0, -1.0}}}},
      {{{std::complex<double>{2.0, -2.0}}}},
  };
  const xt::xtensor<double, 1> kFrequencies{42.0e6};
  const xt::xtensor<int, 2> kNms{{6, 7, 8}};
  const std::array<double, 3> kDirection1{20.0, -10.0, 5.0};
  const std::array<double, 3> kDirection2{3.0, 42.0, -1.0};

  const std::string kFilename = "test_get_model.h5";
  CreateH5File(kFilename, kCoefficients, kFrequencies, kNms);
  const auto element_response =
      std::make_shared<TestResponse>(kFilename, std::nullopt);
  element_response->SetModel(ElementResponseModel::kAartfaacOuter);

  std::shared_ptr<ElementResponse> fixed_direction1 =
      element_response->FixateDirection(kDirection1);
  std::shared_ptr<ElementResponse> fixed_direction2 =
      fixed_direction1->FixateDirection(kDirection2);

  BOOST_TEST(fixed_direction1->GetModel() ==
             ElementResponseModel::kAartfaacOuter);
  BOOST_TEST(fixed_direction2->GetModel() ==
             ElementResponseModel::kAartfaacOuter);

  element_response->SetModel(ElementResponseModel::kOSKARSphericalWave);
  BOOST_TEST(fixed_direction1->GetModel() ==
             ElementResponseModel::kOSKARSphericalWave);
  BOOST_TEST(fixed_direction2->GetModel() ==
             ElementResponseModel::kOSKARSphericalWave);
}

BOOST_AUTO_TEST_SUITE_END()
