#!/bin/bash
./pfDeconv -ref tests/testData/refCountForTesting.csv -alt tests/testData/altCountForTesting.csv -plaf tests/testData/plafForTesting.csv -panel tests/testData/panelForTesting.csv -o tmp
cut -f 1 tmp.hap | tail -n+2 > tmpCHROM
cut -f 2 tmp.hap | tail -n+2 > tmpPOS


awk -F "\"*,\"*" '{print $1}' tests/testData/refCountForTesting.csv  | tail -n+2 > trueCHROM
awk -F "\"*,\"*" '{print $2}' tests/testData/refCountForTesting.csv  | tail -n+2 > truePOS


diff trueCHROM tmpCHROM
if [ $? -ne 0 ]; then
  echo ""
  echo "CHROM unequal"
  exit 1
fi

diff truePOS tmpPOS
if [ $? -ne 0 ]; then
  echo ""
  echo "POS unequal"
  exit 1
fi

rm tmpCHROM tmpPOS trueCHROM truePOS
