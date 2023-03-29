#ifndef LBM_FILE_WRITER_HPP
#define LBM_FILE_WRITER_HPP

#include <fmt/core.h>

#include <Eigen/Core>
#include <cstdio>
#include <filesystem>
#include <fstream>

namespace lbm {

class FileWriter {
 public:
  FileWriter() : output_directory_{}, exists_{false} {}

  /**
   * @brief Construct a new File Writer object
   *
   * @param output_directory
   * @throw std::filesystem::filesystem_error
   */
  FileWriter(const std::filesystem::path& output_directory)
      : output_directory_{output_directory}, exists_{false} {
    namespace fs = std::filesystem;
    if (!fs::exists(output_directory)) {
      fs::create_directories(output_directory);
    }
    exists_ = true;
  }

  void set_output_directory(const std::filesystem::path& output_directory) {
    namespace fs = std::filesystem;
    output_directory_ = output_directory;
    fs::create_directories(output_directory);
    exists_ = true;
  }

  /**
   * @brief Write matrix to file
   *
   * @tparam T
   * @param x Matrix
   * @param filename File name
   */
  template <typename T>
  void write(const Eigen::MatrixBase<T>& x,
             const std::string& filename) const noexcept {
    namespace fs = std::filesystem;
    if (!exists_) {
      fs::create_directories(output_directory_);
    }
    const auto p = output_directory_ / fs::path(filename);
    std::ofstream f(p);
    if (!f) {
      fmt::print(stderr, "Error: could not open a file: {}", p.string());
      std::exit(EXIT_FAILURE);
    }
    f << x << std::endl;
  }

 private:
  std::filesystem::path output_directory_;
  bool exists_;
};

}  // namespace lbm

#endif  // LBM_FILE_WRITER_HPP