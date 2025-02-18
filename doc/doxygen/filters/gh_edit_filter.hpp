// Author: Benedict Armstrong
// Date: February 2025
// This is part of the LehrFEM++ code suite

#pragma once

#include <fstream>
#include <iostream>
#include <optional>
#include <regex>
#include <string>
#include <string_view>

#include "filter.hpp"

/**
 * This Filter replaces all the @gh_edit commands in the source file with a link
 * to the corresponding file on GitHub to allow for easy editing of the
 * documentation.
 */
class GithubEditFilter : public Filter {
 private:
  std::string location_;

  const std::string repo_url =
      "https://github.com/craffael/lehrfempp/blob/master/";

 public:
  explicit GithubEditFilter(const std::string file_path,
                            const std::string project_root);

  ~GithubEditFilter() override = default;

  /**
   * @brief Goes through the source file passed by doxygen (second
   * commandline argument) and replace @gh_edit with the a link to the
   * corresponding file on GitHub.
   */
  std::string filter(const std::string& input) override;
};