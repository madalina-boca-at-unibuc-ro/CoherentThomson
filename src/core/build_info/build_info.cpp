#include "../include/build_info/build_info.hpp"

#include <unistd.h>

#include <fstream>
#include <stdexcept>

namespace Core::BuildInfo {

std::string describe() { return std::string(git_commit) + (git_dirty ? "-dirty" : ""); }

std::string run_host() {
  char name[256] = {};
  if (gethostname(name, sizeof(name) - 1) != 0) return "unknown";
  return name;
}

void require_clean_build() {
  if (!git_available) {
    throw std::runtime_error("require_clean_build is set, but the binary was built without git information");
  }
  if (git_dirty) {
    throw std::runtime_error("require_clean_build is set, but the binary was built from uncommitted changes (" +
                             describe() + "). Commit them and rebuild. Changed files:\n" + git_status);
  }
}

void write_source_patch(const std::string& filepath) {
  if (!git_dirty) return;
  std::ofstream file(filepath, std::ios::binary);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file for source patch export: " + filepath);
  }
  file.write(source_patch.data(), static_cast<std::streamsize>(source_patch.size()));
}

}  // namespace Core::BuildInfo
