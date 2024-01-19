/**
 * @file
 * @brief Defines the LehrFemInfo class
 * @author Raffael Casagrande
 * @date   2020-10-22 09:27:12
 * @copyright MIT License
 */

#ifndef INCG8efea66a20694b749ebcac463c2beaae
#define INCG8efea66a20694b749ebcac463c2beaae

#include <string>

namespace lf::base {

/**
 * @brief Provides extra information about this version of LehrFEM++, in
 * particular licensing information.
 */
class LehrFemInfo {
 public:
  /**
   * @brief Get the git sha1 of this commit of LehrFEM++ (agrees with the sha1
   * of the respective commit on github)
   * @return The commit sha1. Can be the empty string, if we cannot determine
   * it.
   */
  static std::string getVersionSha();

  /**
   * @brief Get the date and time of this commit of LehrFEM++
   * @return The date and time in ISO 8601 form. Can be the empty string, if we
   * could not determine this.
   */
  static std::string getVersionDateTime();

  /**
   * @brief Get the name of a git Tag attached to this commit. (Mostly the empty
   * string)
   * @return The git Tag, if there is one. Otherwise the empty string.
   *
   * @note This function may return an empty string, even if the commit of this
   * version of LehrFEM++ has a tag attached on github. There is no guarantee,
   * that this tag will be returned here, but it might.
   */
  static std::string getVersionTag();

  /**
   * @brief Prints the LehrFEM++ banner to the given stream including version
   * information.
   * @param stream The stream to which the banner should be written.
   */
  static void PrintInfo(std::ostream& stream);

  /**
   * @brief Print the MIT License which applies to the LehrFEM++ code.
   * @param stream The stream to which the license should be output.
   */
  static void PrintLicense(std::ostream& stream);

  /**
   * @brief Print information about third-party libraries used by LehrFEM++ to
   * the given stream.
   * @param stream The stream to which the information should be printed.
   */
  static void Print3rdPartyLicenses(std::ostream& stream);
};

}  // namespace lf::base

#endif  // INCG8efea66a20694b749ebcac463c2beaae
