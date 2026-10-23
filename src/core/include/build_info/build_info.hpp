#pragma once

#include <string>
#include <string_view>

// Provenance of the running binary, captured when it was built (cmake/build_info.cmake), not when it runs: the
// binary may be older than the checkout, or synced from another machine.
namespace Core::BuildInfo {

extern const bool git_available;      // false if git or the repository was unavailable at build time
extern const bool git_dirty;          // uncommitted changes, untracked files included
extern const char* const git_commit;  // full hash of HEAD, or "unknown"
extern const char* const git_commit_date;
extern const char* const git_status;  // `git status --porcelain` lines, empty when clean
extern const char* const compiler;    // compiler id and version
extern const char* const cxx_flags;   // the effective CXX flags of the build type
extern const char* const build_type;
extern const char* const build_host;
extern const std::string_view source_patch;  // the uncommitted changes against git_commit, empty when clean

// "<commit>" or "<commit>-dirty".
std::string describe();

// Hostname of the machine running the binary.
std::string run_host();

// Throws unless the build came from a clean, known commit. Called when the config sets require_clean_build.
void require_clean_build();

// Writes source_patch to `filepath` (`git apply` it on a checkout of git_commit). Does nothing for a clean build.
void write_source_patch(const std::string& filepath);

}  // namespace Core::BuildInfo
