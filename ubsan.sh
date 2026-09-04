#!/bin/sh
# Build ravetools under clang's UndefinedBehaviorSanitizer with -std=gnu++11
# and run the workload in ubsan.R.
#
# Why C++11? `VCGLIB` is written and compiled with the old compilers. Some
# usages were undefined back then, such as left bit shift of `NA_INTEGER`.
# However, `VCGLIB` implements that regardless. This left potential
# vulnerability when the intended undefined behavior is consistent across
# C++ versions. For example, left bit shift of `NA_INTEGER` is well defined
# in C++20, but with totally different behavior. The UBSAN error in C++11
# becomes an implementation error: if we run UBSAN with C++20, there is no way
# we find this bug.
#
# Why C++11 and not the local default: R builds this package with
# -std=gnu++20, where a left shift of a negative value is *well defined*.
# C++11/14/17 all leave it undefined, so only an older-standard build has the
# sensitivity CRAN's clang-UBSAN farm has. That difference is exactly what hid
# the vcglib bit-flag leak fixed in 0.3.1 from every local check.
#
# Mechanism: R's Makeconf sets `CXX = clang++ ... -std=gnu++20` and
# `SHLIB_CXXLD = $(CXX)`, so overriding CXX alone changes the standard *and*
# puts the sanitizer on both the compile and the link line. The flags cannot
# ride on CXXFLAGS because SHLIB_CXXLDFLAGS does not include it.
#
# Usage:  sh ubsan.sh
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="${TMPDIR:-/tmp}/ravetools-ubsan"
SANLIB="$WORK/lib"
mkdir -p "$SANLIB"

# R_MAKEVARS_USER keeps this out of ~/.R/Makevars entirely.
cat > "$WORK/Makevars" <<'EOF'
CXX = clang++ -std=gnu++11 -fsanitize=undefined -fno-sanitize=float-divide-by-zero,vptr -fno-omit-frame-pointer
CXXFLAGS = -g -O1 -Wall
EOF

echo "==> building into $SANLIB"
# --preclean is REQUIRED: there is no header dependency tracking here, so
# stale .o files otherwise survive header-only edits (vcglib is header-only).
R_MAKEVARS_USER="$WORK/Makevars" \
  R CMD INSTALL --preclean --no-multiarch -l "$SANLIB" "$HERE"

echo "==> running workload"
# halt_on_error=0 so one run surfaces many findings; UBSAN already reports each
# source location only once, so this does not flood.
# print_stacktrace is omitted deliberately: llvm-symbolizer is not installed,
# so frames would be raw addresses. The "runtime error: file:line:col" line is
# emitted from the instrumented source itself and needs no symbolizer.
UBSAN_OPTIONS=halt_on_error=0:report_error_type=1 \
RAVETOOLS_SRC="$HERE" \
R_LIBS="$SANLIB" \
  Rscript "$HERE/ubsan.R" 2>&1 | tee "$WORK/ubsan.log"

echo
echo "==================== UBSAN findings ===================="
if grep -E "runtime error:" "$WORK/ubsan.log" > "$WORK/findings.txt" 2>/dev/null; then
  sort -u "$WORK/findings.txt"
  echo
  echo "$(sort -u "$WORK/findings.txt" | wc -l | tr -d ' ') unique finding(s); full log: $WORK/ubsan.log"
  exit 1
else
  echo "none"
  echo "full log: $WORK/ubsan.log"
fi
