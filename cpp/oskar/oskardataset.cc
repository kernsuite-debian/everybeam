// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <iostream>
#include <cassert>
#include <stdexcept>

#include "oskardataset.h"

Dataset::Dataset(H5::H5File& h5_file, const unsigned int freq) {
  // Try to read coefficients for this frequency
  std::string dataset_name = std::to_string((int)(freq / 1e6));

  try {
    // Open dataset
    H5::DataSet dataset = h5_file.openDataSet(dataset_name);

    // Read dataset dimensions
    H5::DataSpace dataspace = dataset.getSpace();
    unsigned int rank = dataspace.getSimpleExtentNdims();
    assert(rank == dataset_rank_);

    // Get dimensions
    std::vector<hsize_t> dims(rank);
    dataspace.getSimpleExtentDims(dims.data(), NULL);
    auto nr_elements = dims[0];
    auto nr_coeffs = dims[1];
    assert(dims[2] == 4);  // tetm*pol

// Coefficient data stored as:
// [nr_elements][nr_coefficients][4],
// with inner dimension:
// (x_te_re, x_te_im), (x_tm_re, x_tm_im),
// (y_te_re, y_te_im), (y_tm_re, y_tm_im)
#ifndef NDEBUG
    std::cout << "nr_elements: " << nr_elements << std::endl;
    std::cout << "nr_coeffs: " << nr_coeffs << std::endl;
#endif

    // Check total number of coefficients to find l_max
    auto l_max_d = sqrt(nr_coeffs + 1) - 1;
    auto l_max = (int)round(l_max_d);
#ifndef NDEBUG
    std::cout << "l_max: " << l_max << std::endl;
#endif

    // Sanity check
    assert(l_max * (l_max + 2) == nr_coeffs);

    // Set members
    nr_elements_ = nr_elements;
    nr_coeffs_ = nr_coeffs;
    l_max_ = l_max;

    // Read coefficients into data vector
    data_.resize(nr_elements * nr_coeffs * 4);
    assert(dims[0] * dims[1] * dims[2] == data_.size());
    H5::DataType data_type = dataset.getDataType();
    assert(data_type.getSize() == sizeof(std::complex<double>));
    dataset.read(data_.data(), data_type, dataspace);
  } catch (H5::FileIException& e) {
    std::stringstream message;
    message << "Could not load dataset for frequency " << dataset_name
            << " Mhz";
    throw std::runtime_error(message.str());
  }
}

size_t Dataset::GetIndex(const unsigned int element) const {
  return element * nr_coeffs_ * 4;
}

std::complex<double>* Dataset::GetAlphaPtr(const unsigned int element) {
  assert(element < GetNrElements());
  size_t index = GetIndex(element);
  return data_.data() + index;
}
