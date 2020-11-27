/**
 * @file
 * @brief Shows how to use loggers to output additional information to the
 * console when refining a mesh.
 * @author Raffael Casagrande
 * @date   2020-10-14 03:40:53
 * @copyright MIT License
 */

// include these headers to load log levels from environment variables/command
// line arguments
#include <spdlog/cfg/argv.h>
#include <spdlog/cfg/env.h>
#include <spdlog/spdlog.h>

#include <spdlog/sinks/basic_file_sink.h>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <iostream>

#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/refinement/refinement.h"

/**
 * @brief Example that demonstrates how spdlog loggers in LehrFEM++ can be used
 * to show additional information.
 */
int main(int argc, char** argv) {
  namespace po = boost::program_options;

  // configure command line arguments:
  po::options_description desc(
      R"foo(This example demonstrates how the loggers of LehrFEM++ can be used to output at different log levels and to different output media.
You can set the log level of every logger via the environment variable SPDLOG_LEVEL or the commandline argument SPDLOG_LEVEL.
Example:
  // set log level of all loggers to debug, set the one of MeshHierarchy to trace:
  export SPDLOG_LEVEL="debug,lf::refinement::MeshHierarchy::Logger=trace"
  ./examples.logger.mesh_hierarchy_demo

  // disable all loggers except for the MeshHierarchy one:
  ./examples.logger.mesh_hierarchy_demo SPDLOG_LEVEL=Off,lf::refinement::MeshHierarchy::Logger=info

For every logger, you can set the level to one of [trace, debug, info, warning, error, fatal, critical, off].

Other available command line options)foo");

  // clang-format off
  desc.add_options()
    ("help,h", "show this help message")
    ("list-loggers,l", "list all loggers that are registered with their log level")
    ("to-file", po::value<std::vector<std::string>>(), "Value of the form '<logger>=<filename>'. Redirects the output of 'logger' to a file name 'filename' instead of standard output.")
  ;

  // clang-format on
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);

  // print help message if needed
  if (vm.count("help") > 0) {
    std::cout << desc << std::endl;
    return 1;
  }

  // load log levels of loggers from environment variable SPDLOG_LEVEL:
  spdlog::cfg::load_env_levels();

  // load log levels from command line arguments:
  spdlog::cfg::load_argv_levels(argc, argv);

  if (vm.count("list-loggers") > 0) {
    // go through all registered loggers and print them with their log level:
    std::cout << "The following loggers are known (with their log level):"
              << std::endl;

    spdlog::apply_all([&](const std::shared_ptr<spdlog::logger>& logger) {
      if (logger->name().empty()) {
        return;  // ignore default logger
      }
      std::cout << logger->name() << ": " << logger->level() << std::endl;
    });

    return 1;
  }

  if (vm.count("to-file") > 0) {
    // redirect loggers to the given file:
    for (const auto& s : vm["to-file"].as<std::vector<std::string>>()) {
      std::vector<std::string> splitted;
      boost::split(splitted, s, boost::is_any_of("="));
      if (splitted.size() != 2) {
        std::cout << "Error, the argument '" << s
                  << "' (given to --to-file) has invalid format" << std::endl;
        return 1;
      }

      // Retrieve logger from registry
      auto logger = spdlog::get(splitted[0]);
      if (logger == nullptr) {
        std::cout << "Error: Could not find logger with name " << splitted[0]
                  << std::endl;
        return 1;
      }

      // replace current sink with a file sink:
      LF_ASSERT_MSG(logger->sinks().size() == 1,
                    "Something is wrong, logger "
                        << logger->name() << " doesn't have exactly 1 sink.");
      logger->sinks()[0] = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
          splitted[1], true);
    }
  }

  // generate a test mesh and refine it once:
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(7);
  lf::refinement::MeshHierarchy mh(
      mesh, std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
  mh.RefineRegular();
}
