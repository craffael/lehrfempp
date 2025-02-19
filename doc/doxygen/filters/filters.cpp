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

// Read file into string
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

  // Define and parse the program options
  std::string aux_file_name = "NPDEFLrefs.aux";
  std::string file_name;
  std::string project_root = "";

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce help message")(
      "aux,a",
      po::value<std::string>(&aux_file_name)->default_value("NPDEFLrefs.aux"),
      "auxiliary file")(
      "project-root,r",
      po::value<std::string>(&project_root)->default_value(""),
      "project root directory")(
      "input-file", po::value<std::string>(&file_name), "input file");

  po::positional_options_description p;
  p.add("input-file", -1);

  po::variables_map vm;
  try {
    po::store(
        po::command_line_parser(argc, argv).options(desc).positional(p).run(),
        vm);
    po::notify(vm);
  } catch (const po::error& ex) {
    std::cerr << ex.what() << std::endl;
    return -1;
  }

  if (vm.count("help")) {
    std::cout << "Usage: " << argv[0] << " [options] <source file>\n";
    std::cout << desc << "\n";
    return 0;
  }

  if (!vm.count("input-file")) {
    std::cerr << "Usage: " << argv[0] << " [options] <source file>\n";
    return -1;
  }

  auto file_string_ = readFile(file_name);

  // Read file into memory
  if (!file_string_.has_value()) {
    std::cerr << "filters.cpp : Could not open " << file_name << std::endl
              << std::flush;
    std::exit(1);
  }

  // Filters
  std::unique_ptr<Filter> filters[] = {
      std::make_unique<GithubEditFilter>(file_name, project_root),
      std::make_unique<LectureDocRefFilter>(file_name, aux_file_name)};

  std::string result = file_string_.value();

  for (auto& filter : filters) {
    result = filter->filter(result);
  }

  std::cout << result;
  return 0;
}