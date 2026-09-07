#!/usr/bin/env bash

set -euo pipefail

usage() {
  printf 'Usage: %s\n' "${0##*/}"
}

case "${1-}" in
  "")
    ;;
  -h|--help)
    usage
    exit 0
    ;;
  *)
    printf '%s: unexpected argument: %s\n' "${0##*/}" "$1" >&2
    usage >&2
    exit 2
    ;;
esac

script_dir=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH='' cd -- "$script_dir/.." && pwd)
cd "$repo_root"

required_commands=(Rscript g++ clang++ clang-format)
for required_command in "${required_commands[@]}"; do
  if ! command -v "$required_command" >/dev/null 2>&1; then
    printf '%s: required command not found: %s\n' \
      "${0##*/}" "$required_command" >&2
    exit 1
  fi
done

required_formatter_version=21.1.8
formatter_version_output=$(clang-format --version)
if [[ $formatter_version_output =~ clang-format\ version\ ([0-9]+\.[0-9]+\.[0-9]+) ]]; then
  formatter_version=${BASH_REMATCH[1]}
else
  printf '%s: unable to determine clang-format version from: %s\n' \
    "${0##*/}" "$formatter_version_output" >&2
  exit 1
fi
if [[ $formatter_version != "$required_formatter_version" ]]; then
  printf '%s: clang-format %s is required, found %s\n' \
    "${0##*/}" "$required_formatter_version" "$formatter_version" >&2
  exit 1
fi

r_include=$(Rscript --vanilla -e 'cat(R.home("include"))')
if [[ ! -d $r_include ]]; then
  printf '%s: R include directory not found: %s\n' \
    "${0##*/}" "$r_include" >&2
  exit 1
fi

mapfile -t linking_to_packages < <(
  Rscript --vanilla -e '
    description <- read.dcf("DESCRIPTION")
    packages <- trimws(strsplit(description[1L, "LinkingTo"], ",")[[1L]])
    cat(packages, sep = "\n")
  '
)
if ((${#linking_to_packages[@]} == 0)); then
  printf '%s: DESCRIPTION does not declare any LinkingTo packages\n' \
    "${0##*/}" >&2
  exit 1
fi

dependency_include_flags=()
for linking_to_package in "${linking_to_packages[@]}"; do
  package_include=$(
    Rscript --vanilla -e '
      package <- commandArgs(trailingOnly = TRUE)[[1L]]
      cat(system.file("include", package = package))
    ' "$linking_to_package"
  )
  if [[ ! -d $package_include ]]; then
    printf '%s: include directory unavailable for LinkingTo package %s\n' \
      "${0##*/}" "$linking_to_package" >&2
    exit 1
  fi
  dependency_include_flags+=(-isystem "$package_include")
done

maintained_sources=(
  src/mutual-neighbor.cpp
  src/neighbor-overlap.cpp
  src/random-dist.cpp
  src/rnx.cpp
  src/triplet-eval.cpp
)

format_sources=(
  "${maintained_sources[@]}"
  src/distance.h
  src/native-validation.h
)

compiler_flags=(
  -std=gnu++17
  -pthread
  -fPIC
  -DRCPP_NO_RTTI
  -DSTRICT_R_HEADERS
  -DRCPP_NO_MODULES
  -Isrc
  -isystem "$r_include"
  "${dependency_include_flags[@]}"
  -Wall
  -Wextra
  -Wpedantic
  -Wformat=2
  -Wnull-dereference
  -Werror
  -fsyntax-only
)

compilers=(g++ clang++)
for compiler in "${compilers[@]}"; do
  for source_file in "${maintained_sources[@]}"; do
    "$compiler" "${compiler_flags[@]}" "$source_file"
  done
done

clang-format --dry-run --Werror "${format_sources[@]}"

printf 'native-quality: PASS (GCC, Clang, clang-format)\n'
