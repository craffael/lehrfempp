#include "gh_edit_filter.hpp"

GithubEditFilter::GithubEditFilter(std::string file_name)
    : file_name_(std::move(file_name)) {}

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

    result += "<a href='" + repo_url + file_name_ +
              "' title='Edit Page on GitHub'>Edit on GH</a>";

    begin += match.prefix().length() + match[0].length();
  }
  result.append(begin, input.cend());

  return result;
}
