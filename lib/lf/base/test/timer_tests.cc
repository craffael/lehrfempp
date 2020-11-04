/**
 * @file
 * @brief Test the implementation of lf::base::Timer and lf::base::AutoTimer
 * @author Raffael Casagrande
 * @date   2020-10-30 02:48:19
 * @copyright MIT License
 */

#include <thread>

#include <gtest/gtest.h>
#include <lf/base/base.h>
#include <spdlog/spdlog.h>

using namespace std::chrono_literals;

namespace lf::base::test {

TEST(lf_base_timer, TimerSingleThread) {
  Timer t;
  // we just sleep -> only wall time progresses, cpu + user time do not progress
  std::this_thread::sleep_for(std::chrono::seconds(1));
  t.Stop();

  EXPECT_TRUE(t.IsStopped());
  auto elapsed = t.Elapsed();
  std::cout << t.Format() << std::endl;
  EXPECT_TRUE(abs(elapsed.wall - 1s) < 0.2s);
  EXPECT_EQ(elapsed.system, 0ns);
  EXPECT_EQ(elapsed.user, 0ns);

  t.Resume();
  EXPECT_FALSE(t.IsStopped());
  // do a real operation that uses mostly system time (allocation + start a
  // thread):
  std::vector<int> x;
  std::thread t3([&]() { x.resize(40000000); });
  std::vector<std::thread> temp;
  for (int i = 0; i < 200; ++i) {
    temp.emplace_back([]() {});
  }
  t3.join();
  for (int i = 0; i < 200; ++i) {
    temp[i].join();
  }
  t.Stop();
  std::cout << t.Format() << std::endl;
  auto elapsed2 = t.Elapsed();
  EXPECT_GT(elapsed2.wall, elapsed.wall);
  EXPECT_GT(elapsed2.system, elapsed.system);
  EXPECT_GE(elapsed2.user, elapsed.user);

  // The following in principle makes sense, but is often not true on
  // travis/github actions:
  // EXPECT_GT(elapsed2.system, elapsed2.user);

  t.Start();
  EXPECT_LT(t.Elapsed().wall, 1ms);
  EXPECT_LT(t.Elapsed().system, 1ms);
  EXPECT_LT(t.Elapsed().user, 1ms);
  for (int i = 2; i < x.size(); ++i) {
    x[i] = x[i - 1] + x[i - 2];
  }
  t.Stop();
  std::cout << t.Format() << std::endl;
  auto elapsed3 = t.Elapsed();
  EXPECT_LT(elapsed3.wall, 3s);
  EXPECT_LT(elapsed3.system, 0.1 * elapsed3.user);
  EXPECT_GT(elapsed3.user, 1ms);
  EXPECT_LT(elapsed3.user, 3s);
}

TEST(lf_base_timer, TimerMultiThread) {
  Timer t(false);
  EXPECT_TRUE(t.IsStopped());
  auto elapsed = t.Elapsed();
  EXPECT_EQ(elapsed.system, 0ns);
  EXPECT_EQ(elapsed.user, 0ns);
  EXPECT_EQ(elapsed.wall, 0ns);

  t.Start();
  // do some multithread work -> system + user > wall
  auto f = []() {
    std::vector<int> x(20000000);
    for (int i = 3; i < x.size(); ++i) {
      x[i] = x[i - 1] + x[i - 2] + x[i - 3];
    }
  };
  std::thread t1(f), t2(f);
  t1.join();
  t2.join();
  t.Stop();
  std::cout << t.Format() << std::endl;
  elapsed = t.Elapsed();
  EXPECT_GT(elapsed.system + elapsed.user, 1.1 * elapsed.wall);
}

TEST(lf_base_timer, AutoTimerSS) {
  std::stringstream ss;
  {
    AutoTimer at{ss};
    auto elapsed = at.Elapsed();
    EXPECT_LT(elapsed.wall, 1ms);
    EXPECT_LT(elapsed.system, 1ms);
    EXPECT_LT(elapsed.user, 1ms);

    std::vector<int> x(20000000);
    EXPECT_EQ(&ss, &at.ostream());
    EXPECT_THROW({ auto x = at.logger(); }, std::bad_variant_access);
  }
  EXPECT_GT(ss.str().size(), 4);
}

TEST(lf_base_timer, AutoTimerLogger) {
  {
    AutoTimer at{spdlog::get(""), spdlog::level::debug};
    auto elapsed = at.Elapsed();
    EXPECT_LT(elapsed.wall, 1ms);
    EXPECT_LT(elapsed.system, 1ms);
    EXPECT_LT(elapsed.user, 1ms);

    std::vector<int> x(20000000);
    EXPECT_EQ(spdlog::get(""), at.logger());
    EXPECT_THROW({ auto& x = at.ostream(); }, std::bad_variant_access);
  }
}

}  // namespace lf::base::test
