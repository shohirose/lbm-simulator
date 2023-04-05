#ifndef LBM_FILE_WRITER_HPP
#define LBM_FILE_WRITER_HPP

#include <fmt/core.h>

#include <Eigen/Core>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace lbm {

class FileWriter {
 public:
  /**
   * @brief Construct a new File Writer object
   *
   * @param output_directory
   * @throw std::filesystem::filesystem_error
   */
  FileWriter(const std::filesystem::path& output_directory)
      : output_directory_{output_directory} {
    namespace fs = std::filesystem;
    if (!fs::exists(output_directory)) {
      fs::create_directories(output_directory);
    }
  }

  /**
   * @brief Write matrix to file
   *
   * @param x Matrix
   * @param filename File name
   * @throw std::runtime_error When failed to open a file.
   */
  template <typename T>
  void write(const Eigen::MatrixBase<T>& x, const std::string& filename) const {
    namespace fs = std::filesystem;
    const auto path = output_directory_ / fs::path(filename);
    std::ofstream file(path);
    if (!file) {
      throw std::runtime_error(
          fmt::format("Error: could not open a file: {}", path.string()));
    }
    file << x << std::endl;
  }

 private:
  std::filesystem::path output_directory_;
};

}  // namespace lbm

#endif  // LBM_FILE_WRITER_HPP