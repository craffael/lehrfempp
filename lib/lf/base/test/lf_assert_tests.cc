/**
 * @file
 * @brief Test that LF_ASSERT/verify work as expected.
 * @author Raffael Casagrande
 * @date   2021-01-27 08:24:32
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/base/base.h>

#include <boost/assert.hpp>

namespace lf::base {

void LfAssertFails() { LF_ASSERT_MSG(false, "hello assert"); }

void LfVerifyFails() { LF_VERIFY_MSG(false, "hello verify"); }

void BoostAssertFails() { BOOST_ASSERT(false); }

void BoostAssertMsgFails() { BOOST_ASSERT_MSG(false, "hello boost assert"); }

void BoostVerifyFails() { BOOST_VERIFY(false); }

void BoostVerifyMsgFails() { BOOST_VERIFY_MSG(false, "hello boost verify"); }

void EigenAssertFails() {
  Eigen::MatrixXd x(2, 2);
  x(2, 0) = 1.0;
}

TEST(lf_base, LfAssertDeathTest) {
  // Make sure the message is contained in the output:
  ASSERT_DEBUG_DEATH(LfAssertFails(), "hello assert");

  // Make sure the function name is included in the stacktrace:
  ASSERT_DEBUG_DEATH(LfAssertFails(), "LfAssertFails");
}

TEST(lf_base, LfVerifyDeathTest) {
  // Make sure the message is contained in the output:
  ASSERT_DEATH(LfVerifyFails(), "hello verify");

  // Make sure the function name is included in the stacktrace:
  ASSERT_DEATH(LfVerifyFails(), "");
}

TEST(lf_base, BoostAssertDeathTest) {
  ASSERT_DEBUG_DEATH(BoostAssertFails(), "BoostAssertFails");
  ASSERT_DEBUG_DEATH(BoostAssertMsgFails(), "BoostAssertMsgFails");
  ASSERT_DEBUG_DEATH(BoostAssertMsgFails(), "hello boost assert");
}

TEST(lf_base, BoostVerifyDeathTest) {
  ASSERT_DEBUG_DEATH(BoostVerifyFails(), "BoostVerifyFails");
  ASSERT_DEBUG_DEATH(BoostVerifyFails(), "");
  ASSERT_DEBUG_DEATH(BoostVerifyMsgFails(), "BoostVerifyMsgFails");
  ASSERT_DEBUG_DEATH(BoostVerifyMsgFails(), "hello boost verify");
}

TEST(lf_base, EigenAssertDeathTest) {
  ASSERT_DEBUG_DEATH(EigenAssertFails(), "EigenAssertFails");
}

}  // namespace lf::base
