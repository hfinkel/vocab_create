//          Copyright (C) Hal Finkel 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)
// SPDX-License-Identifier: BSL-1.0

#include <boost/config/warning_disable.hpp>

#define BOOST_SPIRIT_X3_UNICODE
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/home/x3/binary.hpp>
#include <boost/spirit/home/support/iterators/istream_iterator.hpp>
namespace spirit = boost::spirit;
namespace x3 = spirit::x3;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zstd.hpp>
namespace io = boost::iostreams;

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup.hpp>
namespace logging = boost::log;

#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
namespace asio = boost::asio;

#include <boost/algorithm/string.hpp>
namespace algo = boost::algorithm;

#include <boost/container_hash/hash.hpp>

#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
namespace fs = std::filesystem;

namespace {
// A simple implementation of a count-min sketch (CMS). CMSs are a well-known
// data structure for keeping track of frequent items in a stream of data. For
// more information, see: https://florian.github.io/count-min-sketch/ or
// https://en.wikipedia.org/wiki/Count%E2%80%93min_sketch
class cms {
public:
  cms(std::size_t width, std::size_t depth) : width(width), depth(depth) {
    data.resize(width*depth);

    // Compute the integer log2 of the width (there are many ways to do this,
    // see:
    // https://graphics.stanford.edu/~seander/bithacks.html#IntegerLogObvious)
    // This way, for each has we compute, we know how many per-depth indexes we
    // can extract from it.
    bits_per_hash = 0;
    for (auto w = width; w >>= 1;) ++bits_per_hash;
    keys_per_hash = 8*sizeof(std::size_t)/++bits_per_hash;
  }

  void merge(const cms &other) {
    std::transform(data.cbegin(), data.cend(), other.data.cbegin(),
                   data.begin(), std::plus<>{});
  }

  void increment_count(std::string_view key) {
    std::size_t seed;
    std::size_t seed_usages = 0;
    for (std::size_t d = 0; d < depth; ++d) {
      if (!seed_usages) {
        seed = 0;
        boost::hash_combine(seed, d);
        boost::hash_combine(seed, key);
        seed_usages = keys_per_hash - 1;
      } else {
        seed >>= bits_per_hash;
        --seed_usages;
      }

      ++data[d*width + seed % width];
    }
  }

  std::size_t get_min_count(std::string_view key) const {
    std::size_t count = std::numeric_limits<std::size_t>::max();

    std::size_t seed;
    std::size_t seed_usages = 0;
    for (std::size_t d = 0; d < depth; ++d) {
      if (!seed_usages) {
        seed = 0;
        boost::hash_combine(seed, d);
        boost::hash_combine(seed, key);
        seed_usages = keys_per_hash - 1;
      } else {
        seed >>= bits_per_hash;
        --seed_usages;
      }

      count = std::min(count, data[d*width + seed % width]);
    }

    return count;
  }

  void keep_only_top_n(std::size_t n) {
    std::vector<std::size_t> top_counts(depth*n);
    std::partial_sort_copy(data.begin(), data.end(), top_counts.begin(),
                           top_counts.end(), std::greater<std::size_t>());
    std::size_t last_top_count = top_counts.back();

    for (auto &c : data)
      if (c < last_top_count)
        c = 0;
  }

  std::size_t get_all_counts() const {
    return std::reduce(data.begin(), data.end())/depth;
  }

  std::size_t get_width() const {
    return width;
  }

  std::size_t get_depth() const {
    return depth;
  }

protected:
  std::size_t bits_per_hash;
  std::size_t keys_per_hash;
  std::size_t width;
  std::size_t depth;
  std::vector<std::size_t> data;
};

} // anonymous namespace

int main(int argc, char *argv[]) {
  po::options_description description("Usage");
  description.add_options()
    ("help,h", "Display this help message")
    ("filter-buffer-size", po::value<std::streamsize>()->
      default_value(32*1024), "Filter buffer size for decompression")
    ("log-level", po::value<std::string>(), "Output log level")
    ("threads", po::value<std::size_t>(), "Number of threads")
    ("config,c", po::value<std::string>(), "Configuration file")
    ("max-token-length", po::value<std::size_t>()->
      default_value(10), "Maximum token length")
    ("max-tokens", po::value<std::size_t>()->
      default_value(10000), "Maximum number of tokens to generate")
    ("cms-width-factor", po::value<std::size_t>()->
      default_value(100), "Width of the count-min sketch, multiplicative to the maximum number of tokens")
    ("cms-depth", po::value<std::size_t>()->
      default_value(3), "Depth of the count-min sketch")
    ("output-file", po::value<std::string>(), "Output file")
    ("input-directory", po::value<std::vector<std::string>>(), "Input directory");

  try {
    std::vector<std::string> input_directory_paths;
    std::string output_file_path;

    po::positional_options_description pos_opts;
    pos_opts.add("output-file", 1);
    pos_opts.add("input-directory", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(description).positional(pos_opts).run(), vm);

    if (vm.count("config")) {
      std::string config_file = vm["config"].as<std::string>();
      std::ifstream config_file_s(config_file.c_str());
      if (!config_file_s) {
        std::cerr << "Error: Failed to open the configuration file: " << config_file << "\n";
        return 1;
      }

      po::store(po::parse_config_file(config_file_s, description), vm);
    }

    po::notify(vm);

    if (vm.count("help")) {
      std::cout << description;

      std::cout << "\nThe input directory is iterated recursively for input files. All regular files\nare treated as text files, with automated decompression provided for files\nwith extensions:"
#ifdef HAS_BOOST_IOSTREAMS_BZIP2
                   " .bz2 (bzip2)"
#endif // HAS_BOOST_IOSTREAMS_BZIP2

#ifdef HAS_BOOST_IOSTREAMS_GZIP
                   " .gz (gzip)"
#endif // HAS_BOOST_IOSTREAMS_GZIP

#ifdef HAS_BOOST_IOSTREAMS_LZMA
                   " .lzma (lzma)"
#endif // HAS_BOOST_IOSTREAMS_LZMA

#ifdef HAS_BOOST_IOSTREAMS_ZSTD
                   " .zst (zstd)"
#endif // HAS_BOOST_IOSTREAMS_ZSTD

#ifdef HAS_BOOST_IOSTREAMS_ZLIB
                   " .zz (zlib)"
#endif // HAS_BOOST_IOSTREAMS_ZLIB
                    "\n";

      std::cout << "\nNote:\nOptions that can be provided on the command line can"
                   " also be provided in a\nconfiguration file."
                   " The configuration file uses a line-based 'name = value' syntax.\n";
      return 0;
    }

    // The default buffer size is small. See http://hoad.io/symmetric_filter-example/
    std::streamsize filter_buffer_size = vm["filter-buffer-size"].as<std::streamsize>();

    if (vm.count("input-directory"))
      input_directory_paths = vm["input-directory"].as<std::vector<std::string>>();

    if (input_directory_paths.empty()) {
      std::cerr << "Error: No input directories provided" << "\n";
      std::cerr << description;
      return 1;
    }

    if (vm.count("output-file"))
      output_file_path = vm["output-file"].as<std::string>();

    if (output_file_path.empty()) {
      std::cerr << "Error: No output file provided" << "\n";
      std::cerr << description;
      return 1;
    }

    //  Open the file early so that, if there is a problem, the user knows sooner rather than later.
    std::ofstream vocab_out(output_file_path, std::ios::binary);
    if (!vocab_out) {
      std::cerr << "Error: Unable to open the output vocabulary file: " << output_file_path << "\n";
      return 1;
    }

    static const std::string COMMON_FMT("[%TimeStamp%][%Severity%]: %Message%");
    logging::register_simple_formatter_factory<logging::trivial::severity_level, char>("Severity");

    logging::add_console_log(
      std::cout,
      logging::keywords::format = COMMON_FMT,
      logging::keywords::auto_flush = true
    );

    logging::add_common_attributes();

    logging::trivial::severity_level severity = logging::trivial::info;
    if (vm.count("log-level")) {
      std::string severity_string = vm["log-level"].as<std::string>();
      if (!logging::trivial::from_string(severity_string.c_str(), severity_string.size(), severity)) {
        std::cerr << "Error: Invalid logging severity level: " << severity_string << "\n";
        return 1;
      }
    }

    logging::core::get()->set_filter(
      logging::trivial::severity >= severity
    );

    BOOST_LOG_TRIVIAL(info) << "Log severity level: " << severity;

    std::unique_ptr<asio::thread_pool> pool(
      vm.count("threads") ?
        new asio::thread_pool(vm["threads"].as<std::size_t>()) :
        new asio::thread_pool()
    );

    std::size_t max_token_length = vm["max-token-length"].as<std::size_t>();
    BOOST_LOG_TRIVIAL(debug) << "Maximum token length: " << max_token_length;
    std::size_t max_tokens = vm["max-tokens"].as<std::size_t>();
    BOOST_LOG_TRIVIAL(debug) << "Maximum number of tokens: " << max_tokens;
    std::size_t cms_width_factor = vm["cms-width-factor"].as<std::size_t>();
    std::size_t cms_width = cms_width_factor * max_tokens;
    BOOST_LOG_TRIVIAL(debug) << "Count-min-sketch width: " << cms_width;
    std::size_t cms_depth = vm["cms-depth"].as<std::size_t>();
    BOOST_LOG_TRIVIAL(debug) << "Count-min-sketch depth: " << cms_depth;

    cms substring_stats(cms_width, cms_depth);
    std::mutex substring_stats_m;

    std::unordered_map<std::string, std::size_t, boost::hash<std::string>> vocab;
    std::mutex vocab_m;
    for (int phase = 0; phase < 2; ++phase) {
      std::vector<std::future<void>> futures;

      for (const auto &input_directory_path : input_directory_paths) {
        BOOST_LOG_TRIVIAL(info) << "Processing files under: " << input_directory_path;

        for (const auto &dir_entry : fs::recursive_directory_iterator(input_directory_path)) {
          if (!dir_entry.is_regular_file())
            continue;

          auto process_file = [dir_entry, filter_buffer_size, max_token_length,
                               phase, &substring_stats, &substring_stats_m,
                               &vocab, &vocab_m]() {
            BOOST_LOG_TRIVIAL(info) << "Processing " << dir_entry << " on thread "
              << std::this_thread::get_id();

            io::filtering_istream fs;
            fs.unsetf(std::ios::skipws);

  #ifdef HAS_BOOST_IOSTREAMS_BZIP2
            if (algo::iequals(dir_entry.path().extension().string(), ".bz2")) {
              BOOST_LOG_TRIVIAL(debug) << "Using a bzip2 decompressor for " << dir_entry;
              fs.push(io::bzip2_decompressor{}, filter_buffer_size);
            }
  #endif // HAS_BOOST_IOSTREAMS_BZIP2

  #ifdef HAS_BOOST_IOSTREAMS_GZIP
            if (algo::iequals(dir_entry.path().extension().string(), ".gz")) {
              BOOST_LOG_TRIVIAL(debug) << "Using a gzip decompressor for " << dir_entry;
              fs.push(io::gzip_decompressor{}, filter_buffer_size);
            }
  #endif // HAS_BOOST_IOSTREAMS_GZIP

  #ifdef HAS_BOOST_IOSTREAMS_LZMA
            if (algo::iequals(dir_entry.path().extension().string(), ".lzma")) {
              BOOST_LOG_TRIVIAL(debug) << "Using a lzma decompressor for " << dir_entry;
              fs.push(io::lzma_decompressor{}, filter_buffer_size);
            }
  #endif // HAS_BOOST_IOSTREAMS_LZMA

  #ifdef HAS_BOOST_IOSTREAMS_ZSTD
            if (algo::iequals(dir_entry.path().extension().string(), ".zst")) {
              BOOST_LOG_TRIVIAL(debug) << "Using a zstd decompressor for " << dir_entry;
              fs.push(io::zstd_decompressor{}, filter_buffer_size);
            }
  #endif // HAS_BOOST_IOSTREAMS_ZSTD

  #ifdef HAS_BOOST_IOSTREAMS_ZLIB
            if (algo::iequals(dir_entry.path().extension().string(), ".zz")) {
              BOOST_LOG_TRIVIAL(debug) << "Using a zlib decompressor for " << dir_entry;
              fs.push(io::zlib_decompressor{}, filter_buffer_size);
            }
  #endif // HAS_BOOST_IOSTREAMS_ZLIB

            fs.push(io::file_descriptor_source(dir_entry.path()));

            cms file_substring_stats(substring_stats.get_width(),
                                     substring_stats.get_depth());
            std::unordered_map<std::string, std::size_t, boost::hash<std::string>> file_vocab;

            size_t word_count = 0, char_count = 0, ss_count = 0;
            auto process_word = [&](std::string_view word) {
              ++word_count;
              for (size_t substring_length = 1;
                   substring_length <= max_token_length &&
                    substring_length <= word.size();
                   ++substring_length) {
                for (std::size_t start_offset = 0; start_offset < substring_length; ++start_offset) {
                  // We know that we don't have more characters than bytes in the
                  // string, but we don't really know how many Unicode characters
                  // we have. We can parse again to find them.
                  auto sb = word.begin(), se = word.end();
                  x3::parse(sb, se,
                    x3::omit[ (x3::repeat(start_offset)[x3::unicode::char_]) ] >>
                    *( x3::raw[ (x3::repeat(substring_length)[x3::unicode::char_]) ][ (
                        [&](auto &c) { ++ss_count;
                                       auto ir = x3::_attr(c);
                                       std::string_view s(std::string_view(&*ir.begin(),
                                                            std::distance(ir.begin(),
                                                                          ir.end())));
                                       if (phase == 0)
                                         file_substring_stats.increment_count(s);
                                       else if (substring_stats.get_min_count(s) > 0)
                                         ++file_vocab[std::string(s.begin(), s.end())]; }
                      ) ]
                    )
                  );
                }
              }
            };

            auto process_char = [&](std::string_view ch) {
              ++char_count;

              if (phase == 0)
                file_substring_stats.increment_count(ch);
              else if (substring_stats.get_min_count(ch) > 0)
                ++file_vocab[std::string(ch.begin(), ch.end())];
            };

            spirit::istream_iterator fsb(fs), fse;
            x3::parse(
              fsb, fse,
              *(
                (x3::unicode::space)
                | x3::raw[+x3::unicode::alnum][ (
                  [&](auto &c) { auto ir = x3::_attr(c);
                                 process_word(std::string_view(&*ir.begin(),
                                                std::distance(ir.begin(), ir.end()))); }
                  ) ]
                | x3::raw[x3::unicode::char_][ (
                  [&](auto &c) { auto ir = x3::_attr(c);
                                 process_char(std::string_view(&*ir.begin(),
                                                std::distance(ir.begin(), ir.end()))); }
                  ) ]
                | x3::raw[x3::byte_] // Be resilient to other random bytes in the stream.
              )
            );

            if (phase == 0) {
              BOOST_LOG_TRIVIAL(debug) << "Finished processing " << dir_entry << " on thread "
                << std::this_thread::get_id() << " read " << word_count
                << " words, decomposed into " << ss_count << " substrings, and "
                << char_count << " other characters, " << file_substring_stats.get_all_counts()
                << " total count-min-sketch counts";

              std::lock_guard g(substring_stats_m);
              substring_stats.merge(file_substring_stats);
            } else {
              BOOST_LOG_TRIVIAL(debug) << "Finished processing " << dir_entry << " on thread "
                << std::this_thread::get_id() << " recorded " << file_vocab.size()
                << " substrings";

              std::lock_guard g(vocab_m);
              for (const auto &v : file_vocab)
                vocab[v.first] += v.second;
            }
          };

          auto task = std::packaged_task<void()>(process_file);
          futures.push_back(task.get_future());
          asio::post(*pool, std::move(task));
        }
      }

      // Wait for all tasks to be done.
      for (auto &f : futures)
        f.get();

      if (phase == 0) {
        BOOST_LOG_TRIVIAL(debug) << "Finished processing " << futures.size()
                                 << " files with a total of "
                                 << substring_stats.get_all_counts()
                                 << " total count-min-sketch counts";

        substring_stats.keep_only_top_n(max_tokens);
      } else {
        BOOST_LOG_TRIVIAL(debug) << "Finished processing " << futures.size()
                                 << " files with a total of "
                                 << vocab.size()
                                 << " total collected substrings";
      }
    }

    pool->join();

    std::size_t max_substring_count = 0;
    for (const auto &v : vocab)
      max_substring_count = std::max(max_substring_count, v.second);

    std::vector<std::pair<std::string, float>> sorted_vocab;
    for (const auto &v : vocab)
      sorted_vocab.push_back(std::make_pair(v.first,
        v.second / double(max_substring_count)));

    std::sort(sorted_vocab.begin(), sorted_vocab.end(),
      [](const auto &a, const auto &b) { return a.second > b.second; });

    BOOST_LOG_TRIVIAL(debug) << "Final substring count is "
                             << sorted_vocab.size()
                             << " (will truncate to " << max_tokens << ")";

    // Because the sketch is approximate, we may have ended up with more tokens
    // than are requested in the final vocabulary.
    sorted_vocab.resize(std::min(sorted_vocab.size(), max_tokens));

    for (const auto &v : sorted_vocab)
      vocab_out << v.first << "\t" << v.second << "\n";

    return 0;
  } catch(std::exception &e) {
    std::cerr << "Exception: " << e.what() << "\n";
  }

  return 1;
}

