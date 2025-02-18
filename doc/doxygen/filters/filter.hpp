#pragma once

#include <string>

class Filter {
 public:
  /**
   * @brief Construct a new Filter object.
   */
  Filter() = default;

  /**
   * @brief Destroy the Filter object.
   */
  virtual ~Filter() = default;

  /**
   * @brief Filter the input string.
   *
   * @param input The input string to be filtered.
   * @return std::string The filtered string.
   */
  virtual std::string filter(const std::string& input) = 0;

 protected:
  std::string file_name_;
};