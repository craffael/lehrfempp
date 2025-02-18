#include "gh_edit_filter.hpp"

GithubEditFilter::GithubEditFilter(const std::string file_path,
                                   const std::string project_root) {
  // if the path is an absolute path, remove the project root from the path.
  if (file_path.find(project_root) != std::string::npos) {
    location_ = file_path.substr(project_root.size());
  } else {
    location_ = file_path;
  }
  // remove the first '/' if it exists.
  if (location_.front() == '/') {
    location_.erase(0, 1);
  }
}

std::string GithubEditFilter::filter(const std::string& input) {
  // disable sync with C streams -> faster output to command line.
  std::ios::sync_with_stdio(false);

  // construct result with regex.
  std::string result;
  result.reserve(input.size() + 100);
  std::regex lref("@gh_edit", std::regex_constants::optimize);
  std::smatch match;
  auto begin = input.cbegin();
  while (regex_search(begin, input.cend(), match, lref)) {
    result += match.prefix().str();

    result += "<a class='github_edit' href='" + repo_url + location_ +
              "' title='Edit page on GitHub' target='_blank'>" +
              "Edit on GitHub" + "</a>";

    begin += match.prefix().length() + match[0].length();
  }
  result.append(begin, input.cend());

  return result;
}
