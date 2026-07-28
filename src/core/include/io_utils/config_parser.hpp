#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "config_map.hpp"
#include "io_utils.hpp"
namespace Core::IoUtils {

inline ConfigMap read_config_file(const std::string& filepath, const bool print_debug) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open configuration file: " + filepath);
  }

  ConfigMap config_map;
  std::string line;

  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string key;

    // Extract the key (first word)
    if (!(iss >> key) || key[0] == '#') {
      continue;
    }

    // Grab EVERYTHING left on this line into the value string
    std::string remainder;
    if (std::getline(iss, remainder)) {
      size_t comment_pos = remainder.find('#');
      if (comment_pos != std::string::npos) {
        remainder = remainder.substr(0, comment_pos);
      }
      size_t first_non_space = remainder.find_first_not_of(" \t");
      if (first_non_space != std::string::npos) {
        remainder = remainder.substr(first_non_space);
      }

      std::string value = trim_trailing(remainder);

      // 3. Store the full sanitized value string ("10 cycles") into our map
      if (!value.empty()) {
        config_map[key] = value;
        if (print_debug) {
          std::cout << "Config: " << key << " = [" << value << "]\n";
        }
      }
    }
  }

  return config_map;
}

inline void print_config(const ConfigMap& config) {
  std::cout << "---- Configuration Parameters ----\n";
  for (const auto& entry : config) {
    std::cout << entry.first << " = " << entry.second << "\n";
  }
  std::cout << "----------------------------------\n";
}

}  // namespace Core::IoUtils
