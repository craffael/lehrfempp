/**
 * @file
 * @brief Some utility functions to test the lf::brep::occt module
 * @author Raffael Casagrande
 * @date   2020-11-11 02:50:52
 * @copyright MIT License
 */

#ifndef __15c6adac3fc54dd8b4b887262e606ea6
#define __15c6adac3fc54dd8b4b887262e606ea6

#include <filesystem>

#include "lf/brep/occt/occt.h"

namespace lf::brep::occt::test {

inline std::unique_ptr<OcctBrepModel> LoadModel(std::string_view filename) {
  std::filesystem::path here = __FILE__;
  return std::make_unique<OcctBrepModel>(
      (here.remove_filename() / "brep_models" / filename).string());
}

}  // namespace lf::brep::occt::test

#endif  // __15c6adac3fc54dd8b4b887262e606ea6
