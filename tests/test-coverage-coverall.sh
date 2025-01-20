#!/bin/bash

# Generate gcov output
make -mj

# Generate html report
lcov --base-directory . --directory . --zerocounters -q
make check -mj
lcov --base-directory . --directory . -c -o coverage/lcov.info
# --rc lcov_branch_coverage=1 option will turn on branch check
lcov --remove coverage/lcov.info "/usr*" "src/codeCogs/*" "src/vcf/*" "src/lasso/*" "src/export/*" "src/gzstream/*" "tests/unittest/*" -o coverage/lcov.info # remove output for external libraries
