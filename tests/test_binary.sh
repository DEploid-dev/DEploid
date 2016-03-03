#!/bin/bash

function test_pfDeconv {
  echo -n " pfDeconv $@ "
  for i in `seq 1 10`; do
    echo -n "."

    # Test using pfDeconv self-checks
    ./pfDeconv_dbg $@ -seed $i > /dev/null
    if [ $? -ne 0 ]; then
      echo ""
      echo "Executing \"./pfDeconv_dbg $@ -seed $i\" failed."
      echo "Debug Call: make -mj2 pfDeconv_dbg && ./pfDeconv_dbg $@ -seed $i 2>&1 | less"
      exit 1
    fi

    # Test for memory leaks
    valgrind --error-exitcode=1 --leak-check=full -q ./pfDeconv $@ -seed $i > /dev/null
    if [ $? -ne 0 ]; then
      echo ""
      echo "Valgrind check of \"./pfDeconv $@ -seed $i\" failed."
      exit 1
    fi

  done
  echo " done."
}

echo "Test no input"
./pfDeconv_dbg > /dev/null
if [ $? -ne 0 ]; then
  echo ""
  echo "Executing \"./pfDeconv_dbg \" failed."
  echo "Debug Call: make -mj2 pfDeconv_dbg && ./pfDeconv_dbg | less"
  exit 1
fi

# Test for memory leaks
valgrind --error-exitcode=1 --leak-check=full -q ./pfDeconv > /dev/null
if [ $? -ne 0 ]; then
  echo ""
  echo "Valgrind check of \"./pfDeconv \" failed."
  exit 1
fi

echo "Testing examples"
 test_pfDeconv -ref "labStrains/PG0390_first100ref.txt" -alt "labStrains/PG0390_first100alt.txt" -plaf "labStrains/labStrains_first100_PLAF.txt" -panel "labStrains/lab_first100_Panel.txt" -o tmp1 || exit 1
echo ""
