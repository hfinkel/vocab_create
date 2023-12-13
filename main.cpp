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

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <thread>
namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
  po::options_description description("Usage");
  description.add_options()
    ("help,h", "Display this help message")
    ("filter-buffer-size", po::value<std::streamsize>()->
      default_value(32*1024), "Filter buffer size")
    ("log-severity", po::value<std::string>(), "Output log level")
    ("threads", po::value<std::size_t>(), "Number of threads")
    ("config,c", po::value<std::string>(), "Configuration file")
    ("input-directory", po::value<std::vector<std::string>>(), "Input directory");

  try {
    std::vector<std::string> input_directory_paths;

    po::positional_options_description pos_opts;
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

    static const std::string COMMON_FMT("[%TimeStamp%][%Severity%]: %Message%");
    logging::register_simple_formatter_factory<logging::trivial::severity_level, char>("Severity");

    logging::add_console_log(
      std::cout,
      logging::keywords::format = COMMON_FMT,
      logging::keywords::auto_flush = true
    );

    logging::add_common_attributes();

    logging::trivial::severity_level severity = logging::trivial::info;
    if (vm.count("log-severity")) {
      std::string severity_string = vm["log-severity"].as<std::string>();
      if (!logging::trivial::from_string(severity_string.c_str(), severity_string.size(), severity)) {
        std::cerr << "Error: Invalid logging severity level: " << severity_string << "\n";
        return 1;
      }
    }

    logging::core::get()->set_filter(
      logging::trivial::severity >= severity
    );

    BOOST_LOG_TRIVIAL(info) << "Log Severity: " << severity;

    std::unique_ptr<asio::thread_pool> pool(
      vm.count("threads") ?
        new asio::thread_pool(vm["threads"].as<std::size_t>()) :
        new asio::thread_pool()
    );

    for (const auto &input_directory_path : input_directory_paths) {
      BOOST_LOG_TRIVIAL(info) << "Processing files under: " << input_directory_path;

      for (const auto &dir_entry : fs::recursive_directory_iterator(input_directory_path)) {
        if (!dir_entry.is_regular_file())
          continue;

        asio::post(*pool, [dir_entry, filter_buffer_size]() {
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

          size_t word_count = 0, char_count = 0;
          auto process_word = [&](std::string_view word) {
            // std::cout << "WORD: " << word << "\n";
            ++word_count;
          };

          auto process_char = [&](std::string_view ch) {
            // std::cout << "CHAR: " << ch << "\n";
            ++char_count;
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

          std::cout << "Words: " << word_count << "; Chars: " << char_count << "\n";

          BOOST_LOG_TRIVIAL(debug) << "Finished processing " << dir_entry << " on thread "
            << std::this_thread::get_id();
        });
      }
    }

    pool->join();

    return 0;
  } catch(std::exception &e) {
    std::cerr << "Exception: " << e.what() << "\n";
  }

  return 1;
}

