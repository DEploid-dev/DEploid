#!/bin/bash

while read file; do
  cpplint --filter=-build/include_subdir,-build/header_guard ${file} 2>&1 > lint_tmp
  if [ $? -ne 0 ]; then
    echo "Coding style check:"
    echo "cpplint ${file} failed"
    cat lint_tmp
    exit 1
  fi
done < .ci/checkedList

echo "Coding style check: PASS"
