#include "Config.h"
#include "cfg/test/TemporaryPath.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>
#include <string>
#include <vector>

namespace cfg {
namespace {
// 256 sites of Al with a few solutes and one vacancy.
Config MakeAlloy() {
  auto config = GenerateFCC({4, 4, 4}, Element("Al"));
  config.SetAtomElementTypeAtLattice(0, Element("X"));
  for (size_t lattice_id = 5; lattice_id < config.GetNumAtoms(); lattice_id += 7) {
    config.SetAtomElementTypeAtLattice(lattice_id, Element(lattice_id % 2 == 0 ? "Mg" : "Zn"));
  }
  return config;
}

void ExpectSameAtoms(const Config &expected, const Config &actual) {
  ASSERT_EQ(expected.GetNumAtoms(), actual.GetNumAtoms());
  for (size_t atom_id = 0; atom_id < expected.GetNumAtoms(); ++atom_id) {
    ASSERT_EQ(expected.GetElementAtAtomId(atom_id), actual.GetElementAtAtomId(atom_id)) << "atom " << atom_id;
    const auto &expected_position =
        expected.GetLatticeVector()[expected.GetLatticeIdFromAtomId(atom_id)].GetRelativePosition();
    const auto &actual_position =
        actual.GetLatticeVector()[actual.GetLatticeIdFromAtomId(atom_id)].GetRelativePosition();
    for (const auto dimension: All_Dimensions) {
      ASSERT_NEAR(expected_position[dimension], actual_position[dimension], 1e-12) << "atom " << atom_id;
    }
  }
}
}    // namespace

TEST(ConfigIo, PlainAndCompressedFilesRoundTrip) {
  const test::TemporaryDirectory temporary_directory;
  const auto &directory = temporary_directory.GetPath();
  const auto config = MakeAlloy();
  for (const std::string name: {"alloy.cfg", "alloy.cfg.gz"}) {
    const auto filename = (directory / name).string();
    config.WriteConfig(filename);
    ExpectSameAtoms(config, Config::ReadConfig(filename));
  }
}

// The auxiliary values used to be written to the wrong stream: ahead of the header in a plain file, and as plain
// text inside the gzip stream of a compressed one.
TEST(ConfigIo, AuxiliaryColumnsStayWithTheirAtoms) {
  const test::TemporaryDirectory temporary_directory;
  const auto &directory = temporary_directory.GetPath();
  const auto config = MakeAlloy();
  std::vector<double> marker(config.GetNumAtoms());
  for (size_t atom_id = 0; atom_id < marker.size(); ++atom_id) {
    marker[atom_id] = 1000.0 + static_cast<double>(atom_id);
  }
  const auto plain = (directory / "marked.cfg").string();
  config.WriteExtendedConfig(plain, {{"marker", marker}});

  std::ifstream ifs(plain);
  std::string line;
  while (std::getline(ifs, line) && line.rfind("auxiliary[3]", 0) != 0) {}
  ASSERT_TRUE(ifs.good()) << "the header of the auxiliary column is missing";
  for (size_t atom_id = 0; atom_id < marker.size(); ++atom_id) {
    std::getline(ifs, line);    // mass
    std::getline(ifs, line);    // element
    ASSERT_TRUE(static_cast<bool>(std::getline(ifs, line))) << "atom " << atom_id;
    std::istringstream line_stream(line);
    std::vector<double> columns;
    for (double value; line_stream >> value;) {
      columns.push_back(value);
    }
    ASSERT_EQ(columns.size(), 7u) << "atom " << atom_id << ": " << line;
    ASSERT_DOUBLE_EQ(columns.back(), marker[atom_id]);
  }

  const auto compressed = (directory / "marked.cfg.gz").string();
  config.WriteExtendedConfig(compressed, {{"marker", marker}});
  ExpectSameAtoms(config, Config::ReadConfig(compressed));
}

TEST(ConfigIo, TruncatedFileIsRejected) {
  const test::TemporaryDirectory temporary_directory;
  const auto &directory = temporary_directory.GetPath();
  const auto filename = (directory / "alloy.cfg").string();
  MakeAlloy().WriteConfig(filename);
  std::filesystem::resize_file(filename, std::filesystem::file_size(filename) / 2);
  EXPECT_THROW(Config::ReadConfig(filename), std::runtime_error);
}

TEST(ConfigIo, MissingInputAndUnwritableOutputAreReported) {
  const test::TemporaryDirectory temporary_directory;
  const auto &directory = temporary_directory.GetPath();
  EXPECT_THROW(Config::ReadConfig((directory / "nothing.cfg").string()), std::runtime_error);
  const auto config = MakeAlloy();
  EXPECT_THROW(config.WriteConfig((directory / "no_such_directory" / "alloy.cfg").string()), std::runtime_error);
  EXPECT_THROW(config.WriteConfig((directory / "no_such_directory" / "alloy.cfg.gz").string()), std::runtime_error);
  EXPECT_THROW(config.WriteMap((directory / "no_such_directory" / "map0.txt").string()), std::runtime_error);
}
}    // namespace cfg
