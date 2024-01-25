// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "lobeselementresponse.h"

#include <charconv>
#include <cmath>
#include <complex>
#include <filesystem>
#include <map>
#include <optional>
#include <string_view>

#include <aocommon/throwruntimeerror.h>
#include <boost/algorithm/string/predicate.hpp>
#include <H5Cpp.h>

#include "../common/mathutils.h"
#include "../common/sphericalharmonics.h"

#include "config.h"

// There are two main modi for the AARTFAAC telescope, AARTFAAC-6 and
// AARTFAAC-12. To properly use AARTFAAC in LOBEs mode the coefficients of all
// stations need to be available. At the moment of writing only a partial set
// is available. This means only AARTFAAC-6 is tested.
static const std::array<std::string_view, 12> kAartfaacStationNames{
    // Available
    "CS002LBA", "CS003LBA", "CS004LBA", "CS005LBA", "CS006LBA", "CS007LBA",
    "CS001LBA", "CS011LBA", "CS013LBA",
    // Currently unavailable
    "CS017LBA", "CS021LBA", "CS032LBA"};

struct AartfaacStation {
  std::string_view station;
  int element;
};

template <class T>
static T ExtractIntegral(std::string_view string) {
  int value;
  std::from_chars_result result =
      std::from_chars(string.begin(), string.end(), value);
  if (result.ec != std::errc{} || result.ptr != string.end()) {
    aocommon::ThrowRuntimeError("The value '", string,
                                "' can't be converted to a number");
  }
  return value;
}

enum class AartfaacElements { kInner, kOuter };

static std::optional<AartfaacStation> GetAartfaacStation(
    std::string_view station_name, AartfaacElements elements) {
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

std::tuple<std::string, std::optional<std::size_t>>
LOBESElementResponse::SphericalHarmonicsArguments(const std::string& name,
                                                  const Options& options) {
  const std::optional<AartfaacStation> aartfaac_station =
      GetAartfaacStation(name, AartfaacElements::kInner);

  const std::filesystem::path search_path =
      options.coeff_path.empty() ? GetPath("lobes")
                                 : std::filesystem::path(options.coeff_path);
  const std::string_view station_name =
      aartfaac_station ? aartfaac_station->station : name;
  const std::string station_file = "LOBES_" + std::string(station_name) + ".h5";
  const std::filesystem::path coeff_file_path = search_path / station_file;

  std::optional<std::size_t> element_index{std::nullopt};
  if (aartfaac_station) {
    element_index = aartfaac_station->element;
  }

  return {coeff_file_path.string(), element_index};
}

std::shared_ptr<const LOBESElementResponse> LOBESElementResponse::GetInstance(
    const std::string& name, const Options& options) {
  // Using a single LOBESElementResponse object for each name reduces memory
  // usage since the coefficients are only loaded once.
  // Using weak pointers in this map ensures that LOBESElementResponse objects,
  // are deleted when they are no longer used, which saves memory.
  static std::map<std::string, std::weak_ptr<const LOBESElementResponse>>
      name_response_map;
  std::shared_ptr<const LOBESElementResponse> instance;

  auto entry = name_response_map.find(name);
  if (entry == name_response_map.end()) {
    instance = std::make_shared<const LOBESElementResponse>(name, options);
    name_response_map.insert({name, instance});
  } else {
    instance = entry->second.lock();
    if (!instance) {
      instance = std::make_shared<const LOBESElementResponse>(name, options);
      entry->second = instance;
    }
  }
  return instance;
}

}  // namespace everybeam
