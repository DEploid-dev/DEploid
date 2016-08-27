#!/bin/bash

# Generate gcov output
#make -mj

# Generate html report
lcov --base-directory . --directory . --zerocounters -q
make check -mj
lcov --base-directory . --directory . -c -o libbash_test.info
# --rc lcov_branch_coverage=1 option will turn on branch check
lcov --remove libbash_test.info "/usr*" "src/codeCogs/*" "src/random/*" "src/export/*" "src/gzstream/*" "tests/unittest/*" -o libbash_test.info # remove output for external libraries
rm -rf ./testCoverage
genhtml -o ./testCoverage -t "libbash test coverage" --num-spaces 4 libbash_test.info --legend --function-coverage --demangle-cpp
