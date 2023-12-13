#          Copyright (C) Hal Finkel 2023.
# Distributed under the Boost Software License, Version 1.0.
#    (See accompanying file LICENSE_1_0.txt or copy at
#          https://www.boost.org/LICENSE_1_0.txt)
# SPDX-License-Identifier: BSL-1.0

cmake_minimum_required(VERSION 3.19)

project(vocab_create)
set(CMAKE_CXX_STANDARD 17)

find_package(Boost 1.73.0 COMPONENTS container iostreams log program_options REQUIRED)

INCLUDE(CheckCXXSourceCompiles)
foreach(dt bzip2 gzip lzma zlib bzip2 zstd)
  set(rlt "${CMAKE_REQUIRED_LIBRARIES}")
  list(APPEND CMAKE_REQUIRED_LIBRARIES Boost::iostreams)

  string(TOUPPER ${dt} dtu)
  check_cxx_source_compiles("
      #include <boost/iostreams/filter/${dt}.hpp>
      namespace io = boost::iostreams;

      int main() {
        io::${dt}_decompressor c{};
        return 0;
      }" HAS_BOOST_IOSTREAMS_${dtu})

  if (HAS_BOOST_IOSTREAMS_${dtu})
    add_definitions(-DHAS_BOOST_IOSTREAMS_${dtu})
  endif()

  set(CMAKE_REQUIRED_LIBRARIES "${rlt}")
endforeach()

add_executable(vocab_create main.cpp)

get_filename_component(bls "${Boost_LOG_LIBRARY}" LAST_EXT)
if ("${bls}" STREQUAL "${CMAKE_SHARED_LIBRARY_SUFFIX}")
  # Note: This is added to make Boost.Log link correctly by defining BOOST_ALL_DYN_LINK
  # See https://github.com/boostorg/log/issues/46
  target_link_libraries(vocab_create Boost::dynamic_linking)
endif()

target_link_libraries(vocab_create Boost::container Boost::iostreams Boost::log Boost::program_options)

