/**
 * @file
 * @brief Implementation of timer.h
 * @author Raffael Casagrande
 * @date   2020-10-30 01:41:22
 * @copyright MIT License
 */

#include "timer.h"

#include <fmt/format.h>

#include <boost/config.hpp>
#include <chrono>
#include <iostream>

#if defined(_WIN32) || defined(WIN32)
#include <Windows.h>
#else
#include <sys/times.h>
#include <unistd.h>
#endif

namespace /* anonymous */ {

#if defined(_WIN32) || defined(WIN32)
#else
std::int_least64_t tick_factor()  // multiplier to convert ticks
                                  //  to nanoseconds; -1 if unknown
{
  static std::int_least64_t tick_factor = 0;
  if (tick_factor == 0) {
    if ((tick_factor = ::sysconf(_SC_CLK_TCK)) <= 0) {
      tick_factor = -1;
    } else {
      tick_factor = INT64_C(1000000000) / tick_factor;  // compute factor
      if (tick_factor == 0) {
        tick_factor = -1;
      }
    }
  }
  return tick_factor;
}
#endif

void get_cpu_times(lf::base::Timer::cpu_times* const current) {
  std::chrono::nanoseconds x(
      std::chrono::high_resolution_clock::now().time_since_epoch());
  current->wall = std::chrono::nanoseconds(x.count());

#if defined(_WIN32) || defined(WIN32)

  FILETIME creation, exit;
  if (::GetProcessTimes(::GetCurrentProcess(), &creation, &exit,
                        (LPFILETIME)&current->system,
                        (LPFILETIME)&current->user)) {
    current->user *= 100;  // Windows uses 100 nanosecond ticks
    current->system *= 100;
  } else

  {
    current->system = current->user = std::chrono::nanoseconds(-1);
  }
#else
  tms tm;  // NOLINT
  clock_t c = ::times(&tm);
  if (c == static_cast<clock_t>(-1))  // error
  {
    current->system = current->user = std::chrono::nanoseconds(-1);
  } else {
    current->system = std::chrono::nanoseconds(tm.tms_stime + tm.tms_cstime);
    current->user = std::chrono::nanoseconds(tm.tms_utime + tm.tms_cutime);
    int_least64_t factor;
    if ((factor = tick_factor()) != -1) {
      current->user *= factor;
      current->system *= factor;
    } else {
      current->user = current->system = std::chrono::nanoseconds(-1);
    }
  }
#endif
}
}  // namespace

namespace lf::base {
Timer::cpu_times Timer::Elapsed() const noexcept {
  if (IsStopped()) {
    return times_;
  }

  cpu_times current;  // NOLINT
  get_cpu_times(&current);
  current.wall -= times_.wall;
  current.user -= times_.user;
  current.system -= times_.system;

  return current;
}

std::string Timer::Format(std::string_view format) const {
  const double sec = 1000000000.0L;
  double wall_sec = static_cast<double>(times_.wall.count()) / sec;
  auto total = times_.system.count() + times_.user.count();
  double total_sec = static_cast<double>(total) / sec;
  double percent = (total_sec / wall_sec) * 100.0;

  return fmt::format(
      format, fmt::arg("w", wall_sec),
      fmt::arg("u", static_cast<double>(times_.user.count()) / sec),
      fmt::arg("s", static_cast<double>(times_.system.count()) / sec),
      fmt::arg("t", total_sec), fmt::arg("p", percent));
}

void Timer::Start() noexcept {
  is_stopped_ = false;
  get_cpu_times(&times_);
}

void Timer::Stop() noexcept {
  if (IsStopped()) {
    return;
  }
  is_stopped_ = true;

  cpu_times current;  // NOLINT
  get_cpu_times(&current);
  times_.wall = (current.wall - times_.wall);
  times_.user = (current.user - times_.user);
  times_.system = (current.system - times_.system);
}

void Timer::Resume() noexcept {
  if (IsStopped()) {
    cpu_times current(times_);
    Start();
    times_.wall -= current.wall;
    times_.user -= current.user;
    times_.system -= current.system;
  }
}

AutoTimer::AutoTimer(std::string format)
    : format_(std::move(format)), output_(&std::cout) {}

AutoTimer::AutoTimer(std::ostream& stream, std::string format)
    : format_(std::move(format)), output_(&stream) {}

AutoTimer::AutoTimer(std::shared_ptr<spdlog::logger> logger,
                     spdlog::level::level_enum level, std::string format)
    : format_(std::move(format)),
      output_(std::make_pair(std::move(logger), level)) {}

AutoTimer::~AutoTimer() {
  if (!timer_.IsStopped()) {
    timer_.Stop();
    try {
      Report();
    } catch (...) {  // eat any exceptions
    }
  }
}

std::ostream& AutoTimer::ostream() const {
  return *std::get<std::ostream*>(output_);
}

const std::shared_ptr<spdlog::logger>& AutoTimer::logger() const {
  return std::get<std::pair<std::shared_ptr<spdlog::logger>,
                            spdlog::level::level_enum>>(output_)
      .first;
}

const std::string& AutoTimer::FormatString() const { return format_; }

void AutoTimer::Report() {
  auto string = timer_.Format(format_);
  if (auto* ostream = std::get_if<std::ostream*>(&output_)) {
    (**ostream) << string << std::endl;
  } else {
    std::get<1>(output_).first->log(std::get<1>(output_).second, string);
  }
}

Timer::cpu_times AutoTimer::Elapsed() const noexcept {
  return timer_.Elapsed();
}

std::string AutoTimer::Format(std::string_view format) const {
  return timer_.Format(format);
}

}  // namespace lf::base
