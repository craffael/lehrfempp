/**
 * @file
 * @brief Unit tests for the LehrFemInfo class
 * @author Raffael Casagrande
 * @date   2020-10-23 09:29:16
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/base/base.h>

namespace lf::base {

TEST(LehrFemInfo, PrintInfo) { LehrFemInfo::PrintInfo(std::cout); }

TEST(LehrFemInfo, PrintLicense) { LehrFemInfo::PrintLicense(std::cout); }

TEST(LehrFemInfo, Print3rdPartyLicenses) {
  LehrFemInfo::Print3rdPartyLicenses(std::cout);
}

}  // namespace lf::base
