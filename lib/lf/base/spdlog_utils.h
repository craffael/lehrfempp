/**
 * @file
 * @brief SPDLog utilities that are used elsewhere.
 * @author Raffael Casagrande
 * @date   2020-10-08 03:12:05
 * @copyright MIT License
 */

#ifndef INCG6f06a5790b0b46cf94fb3cc3cc0cc2d3
#define INCG6f06a5790b0b46cf94fb3cc3cc0cc2d3

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <spdlog/formatter.h>
#include <spdlog/logger.h>

#include <Eigen/Core>
#include <memory>

namespace lf::base {

/**
 * @brief Create a [spdlog
 * logger](https://github.com/gabime/spdlog/wiki/2.-Creating-loggers), register
 * it in the [spdlog
 * registry](https://github.com/gabime/spdlog/wiki/5.-Logger-registry) and
 * initialize it with LehrFEM++ specific settings.
 * @param name The name of the logger, is usually equal to the fully-scoped
 * variable name, e.g. `lf::mesh::hybrid2d::Mesh::Logger` or
 * `lf::assemble::AssembleMatrixLogger`
 * @return The initialized logger
 *
 * LehrFEM++ uses spdlog loggers to provide the user with
 * additional information about long running operations. These loggers are
 * either associated with a class, in which case they are static, public member
 * functions of this class (e.g. `lf::mesh::hybrid2d::Mesh::Logger()`) or they
 * are associated with a free function in which case they are retrieved through
 * a global free function lying in the same namespace as the free function (e.g.
 * `lf::assemble::AssembleMatrixLogger()`).
 *
 * You should use InitLogger() to initialize these loggers. In this way, we can
 * set some sensible defaults for all LehrFEM++ loggers at one central location.
 *
 * @note So far a call to this method is essentially forwarded to
 * spdlog::stdout_color_mt() and the formatter is set to a LineFeedFormatter.
 */
std::shared_ptr<spdlog::logger> InitLogger(const std::string& name);

// clang-format off
/**
 * @brief A [spdlog
 * formatter](https://github.com/gabime/spdlog/blob/v1.x/include/spdlog/formatter.h)
 * which wraps another formatter and makes sure that if there are new lines
 * (`\n`) in the log message, that the log message is still properly indented.
 *
 * ### Multiline log message with normal pattern formatter:
 * @snippet line_feed_formatter.cc NormalFormatter
 *
 * Output:
 * ```
[2020-10-08 20:45:54.021] [logger_name] [trace] [mesh_factory_tests.cc:77] hello
world
 * ```
 *
 * ### Multiline log message with LineFeedFormatter:
 * @snippet line_feed_formatter.cc LineFeedFormatter
 *
 * Output:
 * ```
[2020-10-08 20:47:15.781] [logger_name] [trace] [mesh_factory_tests.cc:81] hello
                                                                           world
 * ```
 */
// clang-format on
class LineFeedFormatter final : public spdlog::formatter {
 public:
  explicit LineFeedFormatter(
      std::unique_ptr<spdlog::formatter> wrapped_formatter);

  LineFeedFormatter(const LineFeedFormatter&) = delete;
  LineFeedFormatter(LineFeedFormatter&&) = default;
  LineFeedFormatter& operator=(const LineFeedFormatter&) = delete;
  LineFeedFormatter& operator=(LineFeedFormatter&&) = default;
  ~LineFeedFormatter() override = default;

  void format(const spdlog::details::log_msg& msg,
              spdlog::memory_buf_t& dest) override;

  [[nodiscard]] std::unique_ptr<formatter> clone() const override;

 private:
  std::unique_ptr<spdlog::formatter> wrapped_formatter_;
};

namespace internal {
template <class MATRIX, typename = std::enable_if_t<std::is_base_of_v<
                            Eigen::DenseBase<MATRIX>, MATRIX>>>
using enable_if_eigen = MATRIX;
}  // namespace internal

}  // namespace lf::base

/// \cond

/**
 * \brief this is the fmt::formatter which is used to format eigen
 * matrices/arrays
 */
template <class MATRIX>
struct fmt::formatter<lf::base::internal::enable_if_eigen<MATRIX>> {
  constexpr auto parse(const format_parse_context& ctx) {
    const auto* it = ctx.begin();
    const auto* end = ctx.end();

    if (it != end && *it != '}') {
      throw format_error("invalid format");
    }

    return it;
  }

  template <typename FormatContext>
  auto format(const MATRIX& matrix, FormatContext& ctx) {
    std::stringstream ss; // NOLINT(misc-const-correctness)
    ss << matrix.format(clean_fmt);

    auto it = ctx.out();
    auto str = ss.str();
    std::copy(str.begin(), str.end(), it);
    return it;
  }

 private:
  static inline const Eigen::IOFormat clean_fmt =
      Eigen::IOFormat(4, 0, ", ", "\n", "[", "]");
};
/// \endcond

#endif  // INCG6f06a5790b0b46cf94fb3cc3cc0cc2d3
