/**
 * @file
 * @brief Shows usage of lf::base::LineFeedFormatter
 * @author Raffael Casagrande
 * @date   2020-10-08 10:38:42
 * @copyright MIT License
 */

#include <lf/base/base.h>
#include <spdlog/pattern_formatter.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

void NormalFormatter() {
  //! [NormalFormatter]
  std::shared_ptr<spdlog::logger> logger =
      spdlog::stdout_color_mt("logger_name");
  logger->set_level(spdlog::level::trace);
  SPDLOG_LOGGER_TRACE(logger, "hello\nworld");
  //! [NormalFormatter]
}

void LineFeedFormatter() {
  //! [LineFeedFormatter]
  std::shared_ptr<spdlog::logger> logger =
      spdlog::stdout_color_mt("logger_name");
  logger->set_formatter(std::make_unique<lf::base::LineFeedFormatter>(
      std::make_unique<spdlog::pattern_formatter>()));
  logger->set_level(spdlog::level::trace);
  SPDLOG_LOGGER_TRACE(logger, "hello\nworld");
  //! [LineFeedFormatter]
}
