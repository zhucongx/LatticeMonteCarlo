#ifndef LMC_LMC_CFG_TEST_TEMPORARYPATH_HPP_
#define LMC_LMC_CFG_TEST_TEMPORARYPATH_HPP_
// Files and directories for tests. The names contain the process id, so that two builds can run their tests at
// the same time, and everything is removed again unless a test failed.
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <string>
#include <system_error>
#include <unistd.h>

namespace test {
inline std::filesystem::path MakeUniquePath(const std::string &extension) {
  static int counter = 0;
  return std::filesystem::temp_directory_path() /
      ("lmc_test_" + std::to_string(::getpid()) + "_" + std::to_string(counter++) + extension);
}

// What a failed test wrote is kept, to be able to look at it.
inline bool IsWorthKeeping() {
  const auto *test_suite = ::testing::UnitTest::GetInstance()->current_test_suite();
  return test_suite != nullptr && test_suite->Failed();
}

class TemporaryDirectory {
 public:
  TemporaryDirectory()
      : path_(MakeUniquePath("")) {
    std::filesystem::create_directories(path_);
  }
  TemporaryDirectory(const TemporaryDirectory &) = delete;
  TemporaryDirectory &operator=(const TemporaryDirectory &) = delete;
  ~TemporaryDirectory() {
    if (!IsWorthKeeping()) {
      std::error_code error;
      std::filesystem::remove_all(path_, error);
    }
  }
  [[nodiscard]] const std::filesystem::path &GetPath() const {
    return path_;
  }

 private:
  const std::filesystem::path path_;
};

class TemporaryFile {
 public:
  TemporaryFile(const std::string &content, const std::string &extension)
      : path_(MakeUniquePath(extension)) {
    std::ofstream(path_, std::ios::binary) << content;
  }
  TemporaryFile(const TemporaryFile &) = delete;
  TemporaryFile &operator=(const TemporaryFile &) = delete;
  ~TemporaryFile() {
    if (!IsWorthKeeping()) {
      std::error_code error;
      std::filesystem::remove(path_, error);
    }
  }
  [[nodiscard]] std::string GetName() const {
    return path_.string();
  }

 private:
  const std::filesystem::path path_;
};
}    // namespace test
#endif    //LMC_LMC_CFG_TEST_TEMPORARYPATH_HPP_
