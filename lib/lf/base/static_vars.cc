/** @file static_vars.cc */

// Preprocessor flag to prevent double definition of global handler object
// for control variables
#define CTRLVARROOT

#include "static_vars.h"
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

namespace lf::base {

// Handling of global control variables
// Global variable !
StaticVar *ctrl_root = nullptr;

// *******************************************************************************
// GLOBAL FUNCTION: ListCtrlVars = Churns out a complete list of all
// control variables in the global list
// *******************************************************************************
unsigned int ListCtrlVars(std::ostream &out, const StaticVar *ctrl_var_root) {
  // Global list ist default
  if (ctrl_var_root == nullptr) {
    ctrl_var_root = ctrl_root;
  }

  out << "List of flags:" << std::endl;

  unsigned int cnt = 0;
  const StaticVar *root = ctrl_var_root;
  while (root != nullptr) {
    if (!(root->name_).empty()) {
      out << "\tVar '" << root->name_ << "' = " << root->ref_;
      if (!(root->comment_).empty()) {
        out << " [" << root->comment_ << "]" << std::endl;
      } else {
        out << std::endl;
      }
      cnt++;
    }
    root = root->next_;
  }
  return cnt;
}

bool SetCtrlVar(const std::string &varname, int value,
                const StaticVar *ctrl_var_root) {
  // Global list ist default
  if (ctrl_var_root == nullptr) {
    ctrl_var_root = ctrl_root;
  }

  const StaticVar *root = ctrl_var_root;
  while (root != nullptr) {
    if (!(root->name_).empty()) {
      if (root->name_ == varname) {
        root->ref_ = value;
        return true;
      }
    }
    root = root->next_;
  }
  return false;
}

CONTROLDECLARE(read_ctrl_vars_file, "read_ctrl_vars_file");

bool ReadCtrlVarsFile(const std::string &filename,
                      const StaticVar *ctrl_var_root) {
  // Open file for reading
  std::ifstream file(filename, std::ios::in);
  if (file.fail()) {
    return false;
  }

  // Create temporary map for key-value string pairs
  std::map<std::string, std::string> keyval_map;

  std::string keyfield;
  std::string valuefield;
  bool readkey = true;
  do {
    char nextchar = file.get();
    switch (nextchar) {
      case ' ': {
        if (readkey) {
          break;
        }
      }
      case '=': {
        if (readkey) {
          readkey = false;
          valuefield.clear();
          break;
        }
      }
      case '\n':
      case -1: {
        if (read_ctrl_vars_file > 0) {
          std::cout << "Read line: " << keyfield << " -> " << valuefield
                    << std::endl;
        }
        keyval_map[keyfield] = valuefield;
        readkey = true;
        keyfield.clear();
        valuefield.clear();
        break;
      }
      default: {
        if (readkey) {
          keyfield.push_back(nextchar);
        } else {
          valuefield.push_back(nextchar);
        }
      }
    }  // end switch
  }    // end do
  while (!file.eof());

  // Global variable list is default
  if (ctrl_var_root == nullptr) {
    ctrl_var_root = ctrl_root;
  }
  const StaticVar *root = ctrl_var_root;

  while (root != nullptr) {
    if (!(root->name_).empty()) {
      auto keyval_ptr = keyval_map.find(root->name_);
      if (keyval_ptr != keyval_map.end()) {
        if (read_ctrl_vars_file > 0) {
          std::cout << "Setting variable " << root->name_ << std::endl;
        }
        std::istringstream intstr(keyval_ptr->second);
        intstr >> root->ref_;
      } else {
        if (read_ctrl_vars_file > 0) {
          std::cout << "Variable " << root->name_ << " not set" << std::endl;
        }
      }
    }
    root = root->next_;
  }
  return true;
}

}  // end namespace lf::base
