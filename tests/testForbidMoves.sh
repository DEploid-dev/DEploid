#!/bin/bash
function test_foundStatement {
    grep "${testingStateMent}" tmp.dbg > /dev/null
    if [ $? -ne 0 ]; then
      echo ""
      echo ${testingStateMent}
      exit 1
    fi

}

function test_notFoundStatement {
    grep "${testingStateMent}" tmp.dbg > /dev/null
    if [ $? -ne 1 ]; then
      echo ""
      echo ${testingStateMent}
      exit 1
    fi
}

# BaseLine
./pfDeconv_dbg -ref tests/testData/refCountForTesting.csv -alt tests/testData/altCountForTesting.csv -plaf tests/testData/plafForTesting.csv -panel tests/testData/panelForTesting.csv -o tmp > tmp.dbg

testingStateMent="update proportion"
test_foundStatement
testingStateMent="Update Single Hap"
test_foundStatement
testingStateMent="Update Pair Hap"
test_foundStatement

testingStateMent="Update Prop: YES"
test_foundStatement
testingStateMent="Update Single: YES"
test_foundStatement
testingStateMent="Update Pair: YES"
test_foundStatement

testingStateMent="Update Prop: NO"
test_notFoundStatement
testingStateMent="Update Single: NO"
test_notFoundStatement
testingStateMent="Update Pair: NO"
test_notFoundStatement

#Forbid prop
./pfDeconv_dbg -ref tests/testData/refCountForTesting.csv -alt tests/testData/altCountForTesting.csv -plaf tests/testData/plafForTesting.csv -panel tests/testData/panelForTesting.csv -o tmp -forbidUpdateProp > tmp.dbg

testingStateMent="update proportion"
test_notFoundStatement
testingStateMent="Update Single Hap"
test_foundStatement
testingStateMent="Update Pair Hap"
test_foundStatement

testingStateMent="Update Prop: YES"
test_notFoundStatement
testingStateMent="Update Single: YES"
test_foundStatement
testingStateMent="Update Pair: YES"
test_foundStatement

testingStateMent="Update Prop: NO"
test_foundStatement
testingStateMent="Update Single: NO"
test_notFoundStatement
testingStateMent="Update Pair: NO"
test_notFoundStatement

#Forbid Single
./pfDeconv_dbg -ref tests/testData/refCountForTesting.csv -alt tests/testData/altCountForTesting.csv -plaf tests/testData/plafForTesting.csv -panel tests/testData/panelForTesting.csv -o tmp -forbidUpdateSingle > tmp.dbg

testingStateMent="update proportion"
test_foundStatement
testingStateMent="Update Single Hap"
test_notFoundStatement
testingStateMent="Update Pair Hap"
test_foundStatement

testingStateMent="Update Prop: YES"
test_foundStatement
testingStateMent="Update Single: YES"
test_notFoundStatement
testingStateMent="Update Pair: YES"
test_foundStatement

testingStateMent="Update Prop: NO"
test_notFoundStatement
testingStateMent="Update Single: NO"
test_foundStatement
testingStateMent="Update Pair: NO"
test_notFoundStatement

#Forbid Pair
./pfDeconv_dbg -ref tests/testData/refCountForTesting.csv -alt tests/testData/altCountForTesting.csv -plaf tests/testData/plafForTesting.csv -panel tests/testData/panelForTesting.csv -o tmp -forbidUpdatePair > tmp.dbg

testingStateMent="update proportion"
test_foundStatement
testingStateMent="Update Single Hap"
test_foundStatement
testingStateMent="Update Pair Hap"
test_notFoundStatement

testingStateMent="Update Prop: YES"
test_foundStatement
testingStateMent="Update Single: YES"
test_foundStatement
testingStateMent="Update Pair: YES"
test_notFoundStatement

testingStateMent="Update Prop: NO"
test_notFoundStatement
testingStateMent="Update Single: NO"
test_notFoundStatement
testingStateMent="Update Pair: NO"
test_foundStatement

./pfDeconv_dbg -ref tests/testData/refCountForTesting.csv \
-alt tests/testData/altCountForTesting.csv \
-plaf tests/testData/plafForTesting.csv \
-panel tests/testData/panelForTesting.csv \
-o tmp -initialP 0.1 0.2 0.3 0.4 -k 4 -forbidUpdateProp > tmp.dbg

testingStateMent="       0.1	       0.2	       0.3	       0.4"
test_foundStatement
testingStateMent="Initial prob: 0.1 0.2 0.3 0.4"
test_foundStatement

./pfDeconv_dbg -ref tests/testData/refCountForTesting.csv \
-alt tests/testData/altCountForTesting.csv \
-plaf tests/testData/plafForTesting.csv \
-panel tests/testData/panelForTesting.csv \
-o tmp -k 4 -forbidUpdateProp > tmp.dbg

testingStateMent="      0.25	      0.25	      0.25	      0.25"
test_foundStatement
