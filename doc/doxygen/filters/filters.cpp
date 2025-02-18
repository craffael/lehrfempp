// Author: Benedict Armstrong
// Date: February 2025
// This is part of the LehrFEM++ code suite

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "filter.hpp"
#include "gh_edit_filter.hpp"
#include "link_filter.hpp"

namespace po = boost::program_options;

std::optional<std::string> readFile(const std::string& file_name) {
  std::ifstream t(file_name, std::ios::binary);  // open the file
  if (!t.is_open()) {
    std::cerr << "filters.cpp : Could not open " << file_name << std::endl
              << std::flush;
    return {};
  }
  t.seekg(0, std::ios::end);
  std::size_t size = t.tellg();
  std::string result;
  result.resize(size);
  t.seekg(0);
  t.read(&result[0], size);

  return result;
}

int main(int argc, char* argv[]) {
  // takes a file as input and print the content after parsing.
  // the output is cout which can be understood natively by doxygen.
  // Process arguments
  std::string aux_file_name_;
  std::string file_name_;
  switch (argc) {
    case 1: {
      std::cout << "Usage: " << argv[0] << " [aux file] <source file>"
                << std::endl;
      return -1;
    }
    case 2: {
      aux_file_name_ = "NPDEFLrefs.aux";
      break;
    }
    default: {
      aux_file_name_ = argv[1];
      break;
    }
  }

  file_name_ = argv[argc - 1];
  auto file_string_ = readFile(file_name_);

  // Read file into memory

  if (!file_string_.has_value()) {
    std::cerr << "filters.cpp : Could not open " << file_name_ << std::endl
              << std::flush;
    std::exit(1);
  }

  // Filters
  std::unique_ptr<Filter> filters[] = {
      std::make_unique<GithubEditFilter>(file_name_),
      std::make_unique<LectureDocRefFilter>(file_name_, aux_file_name_)};

  std::string result = file_string_.value();

  for (auto& filter : filters) {
    result = filter->filter(result);
  }

  std::cout << result;
  return 0;
}
