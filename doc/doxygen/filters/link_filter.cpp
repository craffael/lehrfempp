#include "link_filter.hpp"

std::optional<std::string> LectureDocRefFilter::ReadFile(
    const std::string& file_name) {
  std::ifstream t(file_name, std::ios::binary);  // open the file
  if (!t.is_open()) {
    return {};
  }
  t.seekg(0, std::ios::end);
  std::size_t aux_size = t.tellg();
  std::string result;
  result.resize(aux_size);
  t.seekg(0);
  t.read(&result[0], aux_size);

  return result;
}

// Constructor
LectureDocRefFilter::LectureDocRefFilter(std::string file, std::string aux_file)
    : file_(std::move(file)), aux_file_(std::move(aux_file)) {
  label_map_ = {{"equation", "Equation"}, {"par", "Paragraph"},
                {"chapter", "Chapter"},   {"sec", "Section"},
                {"section", "Section"},   {"subsection", "Subsection"},
                {"figure", "Figure"},     {"code", "Code"},
                {"remark", "Remark"},     {"subsubsection", "Subsection"},
                {"example", "Example"}};

  link_with_page_num_ = {
      {"Paragraph", true},
      {"Code", true},
      {"Remark", true},
  };
}
void LectureDocRefFilter::LoadAuxTable() {
  namespace qi = boost::spirit::qi;
  namespace karma = boost::spirit::karma;

  // load auxilliary file and compute it's hash
  auto aux_string = ReadFile(aux_file_);
  if (!aux_string.has_value()) {
    std::cerr << "filter.cpp : Could not open " << aux_file_ << std::endl
              << std::flush;
    std::exit(1);
  }
  auto aux_hash = std::hash<std::string>()(aux_string.value());

  // load cache file and retrieve hash
  std::string cache_filename = "filter.cache";
  bool rebuild_cache = false;  // should we save our cache at the end?
  auto cache_string = ReadFile(cache_filename);
  if (cache_string.has_value()) {
    std::size_t cache_hash;
    std::size_t num_entries;
    auto begin = cache_string.value().cbegin();
    auto end = cache_string.value().cend();
    if (!qi::phrase_parse(begin, end, qi::ulong_long >> qi::ulong_long, qi::eol,
                          cache_hash, num_entries)) {
      std::cerr << "filter.cpp : Cannot get hash from " << cache_filename
                << " -> rebuilding cache" << std::endl;
      rebuild_cache = true;
    } else {
      if (cache_hash != aux_hash) {
        // if hash stored in cache file doesn't agree with our computed hash
        // of the aux file, rebuild cache
        rebuild_cache = true;
        std::cerr << "filter.cpp : Hashes don't agree -> Rebuild Cache"
                  << std::endl;
      } else {
        // load the rest of the cache
        aux_table_.reserve(num_entries);
        if (!qi::phrase_parse(begin, end,
                              *(qi::lexeme[+(qi::char_ - qi::eol)] >>
                                qi::lexeme[+(qi::char_ - qi::eol)]),
                              qi::eol, aux_table_)) {
          std::cerr
              << "filter.cpp : Cannot load cache entries -> rebuilding cache";
          rebuild_cache = true;
        }
      }
    }
  } else {
    rebuild_cache = true;
  }

  if (rebuild_cache) {
    // build cache:
    auto begin = aux_string.value().cbegin();

    std::smatch base_match;
    const std::regex pattern(R"(\\newlabel\{([^\@}]*)(?:@cref)?\}\{(.*)\}\n)",
                             std::regex::optimize);
    std::unordered_map<std::string, Match> match_map;

    // Find all lines that at labels
    while (
        regex_search(begin, aux_string.value().cend(), base_match, pattern)) {
      if (base_match.ready() && !base_match.empty()) {
        std::string label = base_match.str(1);
        Match m;

        if (match_map.count(label) > 0) {
          m = match_map[label];
        }

        // Check if label matched the first pattern
        std::smatch sub_match;
        std::string line = base_match.str(0);
        // remove leading and trailing whitespaces
        boost::algorithm::trim(line);

        if (std::regex_match(line, sub_match, label_pattern_1_)) {
          m.Label = label_map_.count(sub_match.str(2)) > 0
                        ? label_map_[sub_match.str(2)]
                        : sub_match.str(2);
          m.number = sub_match.str(3);
          m.page_number = sub_match.str(4);
          match_map[label] = m;
          // std::cerr << "Match1: " << label << ": " << m.Label << " " <<
          // m.number
          //           << " " << m.page_number << std::endl;
        } else if (std::regex_match(line, sub_match, label_pattern_2_)) {
          m.number = sub_match.str(2);
          m.page_number = sub_match.str(3);
          m.title = sub_match.str(4);
          m.reference_name = sub_match.str(5);
          match_map[label] = m;
          // std::cerr << "Match2: " << label << ": " << m.number << " "
          //           << m.page_number << " " << m.title << std::endl;
        } else {
          std::cerr << "filter.cpp : Cannot parse '" << line << "'"
                    << std::endl;
        }

        begin += base_match.prefix().length() + base_match[0].length();
      }
    }

    // Create aux table from match_map
    for (const auto& [label, match] : match_map) {
      if (match.Label.empty() || label.empty()) {
        std::cerr << "filter.cpp : Label empty for " << match.reference_name
                  << std::endl;
        continue;
      }
      std::string s = match.Label;
      s += ";" + (match.number.empty() ? "" : match.number);
      s += ";" + (match.page_number.empty() ? "" : match.page_number);
      s += ";" + (match.title.empty() ? "" : match.title);
      s += ";" + (match.reference_name.empty() ? "" : match.reference_name);

      aux_table_.emplace_back(label, s);
    }

    std::sort(aux_table_.begin(), aux_table_.end(),
              [](auto& a, auto& b) { return a.first < b.first; });

    // save cache
    std::ofstream cache_out(cache_filename, std::ios_base::out |
                                                std::ios_base::binary |
                                                std::ios_base::trunc);
    karma::ostream_iterator<char> outit(cache_out);

    if (!karma::generate_delimited(
            outit,
            karma::ulong_long << karma::ulong_long
                              << *(karma::string << karma::string),
            karma::eol, aux_hash, aux_table_.size(), aux_table_)) {
      std::cerr << "filter.cpp : Error while writing to cache!" << std::endl;
      std::exit(1);
    }
  }
}

std::string LectureDocRefFilter::filter(const std::string& input) {
  // disable sync with C streams -> faster output to command line.
  std::ios::sync_with_stdio(false);

  // construct result with regex.
  std::string result;
  result.reserve(input.size() + 100);
  std::regex lref("@lref(_link)?\\{(.*?)\\}", std::regex_constants::optimize);
  std::smatch match;
  auto begin = input.cbegin();
  while (regex_search(begin, input.cend(), match, lref)) {
    result += match.prefix().str();
    if (aux_table_.empty()) {
      // build aux lookup table only if we need it.
      LoadAuxTable();
    }

    auto search_result =
        std::lower_bound(aux_table_.begin(), aux_table_.end(),
                         std::pair<std::string, std::string>(match[2], ""),
                         [](const std::pair<std::string, std::string>& a,
                            const std::pair<std::string, std::string>& b) {
                           return a.first < b.first;
                         });
    if (search_result == aux_table_.end()) {
      std::cerr << "filter.cpp : Cannot find a label for " << match[2]
                << std::endl;
      result += match[0];
    } else {
      //  Split search_result->second into two parts: Label and link (by #)
      std::string sample = search_result->second;
      std::vector<std::string> strs;
      boost::split(strs, sample, boost::is_any_of(";"));

      // std::cerr << "Match: " << match[0] << "(" << search_result->first <<
      // ","
      //           << sample << ")" << " -> " << strs.size() << ";";
      // for (auto& str : strs) {
      //   std::cerr << str;
      // }
      // std::cerr << strs.size() << std::endl;

      std::string text = strs[0] + " " + strs[1];
      if (match[1].str() == "_link") {
        std::string url = lecture_doc_url_ + "#" + strs[4];
        if (link_with_page_num_.count(strs[0]) > 0 &&
            link_with_page_num_[strs[0]]) {
          url = lecture_doc_url_ + "#page=" + strs[2];
        }

        result += "<a href='" + url + "' title='" + strs[3] +
                  "'>Lecture Document " + text + "</a>";
      } else {
        result += text;
      }
    }
    begin += match.prefix().length() + match[0].length();
  }
  result.append(begin, input.cend());

  return result;
}

// int main(int argc, char* argv[]) {
//   // takes a file as input and print the content after parsing.
//   // the output is cout which can be understood natively by doxygen.
//   // Process arguments
//   std::string aux_file;
//   std::string file;
//   switch (argc) {
//     case 1: {
//       std::cout << "Usage: " << argv[0] << " [aux file] <source file>"
//                 << std::endl;
//       return -1;
//     }
//     case 2: {
//       aux_file = "NPDEFLrefs.aux";
//       file = argv[1];
//       break;
//     }
//     default: {
//       aux_file = argv[1];
//       file = argv[2];
//       break;
//     }
//   }

//   // Perform filtering
//   LectureDocRefFilter f(file, aux_file);
//   std::cout << f.Filter();

//   return 0;
// }
