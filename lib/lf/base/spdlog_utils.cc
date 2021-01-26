/**
 * @file
 * @brief Implementations of spdlog_utils.h
 * @author Raffael Casagrande
 * @date   2020-10-08 03:46:12
 * @copyright MIT License
 */

#include "spdlog_utils.h"

#include <spdlog/pattern_formatter.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace lf::base {

std::shared_ptr<spdlog::logger> InitLogger(const std::string& name) {
  auto result = spdlog::stdout_color_mt(name);
  result->set_formatter(std::make_unique<LineFeedFormatter>(
      std::make_unique<spdlog::pattern_formatter>()));
  return result;
}

LineFeedFormatter::LineFeedFormatter(
    std::unique_ptr<spdlog::formatter> wrapped_formatter)
    : wrapped_formatter_(std::move(wrapped_formatter)) {}

void LineFeedFormatter::format(const spdlog::details::log_msg& msg,
                               spdlog::memory_buf_t& dest) {
  int offset = -2;  // -2: first line, // -1 second line
  const auto* begin = msg.payload.begin();

  for (const auto* it = msg.payload.begin();; ++it) {
    if (it == msg.payload.end() || *it == '\n') {
      if (offset == -2) {
        ++offset;
      } else if (offset == -1) {
        // the log message has more than one line
        // => determine the offset
        spdlog::memory_buf_t temp_dest;
        spdlog::details::log_msg empty_msg(msg.time, msg.source,
                                           msg.logger_name, msg.level,
                                           spdlog::string_view_t());
        wrapped_formatter_->format(empty_msg, temp_dest);
        offset = static_cast<int>(temp_dest.size() - 1);
      }

      spdlog::details::log_msg part(msg.time, msg.source, msg.logger_name,
                                    msg.level,
                                    spdlog::string_view_t(begin, it - begin));

      auto old_end = dest.size();
      wrapped_formatter_->format(part, dest);

      if (offset > 0) {
        std::memset(dest.data() + old_end, ' ', offset);
      }

      if (it == msg.payload.end()) {
        break;
      }
      begin = it + 1;
    }
  }
}

std::unique_ptr<spdlog::formatter> LineFeedFormatter::clone() const {
  return std::make_unique<LineFeedFormatter>(wrapped_formatter_->clone());
}

}  // namespace lf::base
