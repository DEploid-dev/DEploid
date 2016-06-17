#!/bin/bash
./pfDeconv -ref tests/testData/refCountForTesting.csv -alt tests/testData/altCountForTesting.csv -plaf tests/testData/plafForTesting.csv -panel tests/testData/panelForTesting.csv -o tmp -exportPostProb -k 2
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

cut -f 1 tmp.single1 | tail -n+2 > tmpCHROM1
cut -f 2 tmp.single1 | tail -n+2 > tmpPOS1

diff trueCHROM tmpCHROM1
if [ $? -ne 0 ]; then
  echo ""
  echo "CHROM unequal"
  exit 1
fi

diff truePOS tmpPOS1
if [ $? -ne 0 ]; then
  echo ""
  echo "POS unequal"
  exit 1
fi

cut -f 1 tmp.single2 | tail -n+2 > tmpCHROM2
cut -f 2 tmp.single2 | tail -n+2 > tmpPOS2

diff trueCHROM tmpCHROM2
if [ $? -ne 0 ]; then
  echo ""
  echo "CHROM unequal"
  exit 1
fi

diff truePOS tmpPOS2
if [ $? -ne 0 ]; then
  echo ""
  echo "POS unequal"
  exit 1
fi


cut -f 1 tmp.pair | tail -n+2 > tmpCHROM3
cut -f 2 tmp.pair | tail -n+2 > tmpPOS3

diff trueCHROM tmpCHROM3
if [ $? -ne 0 ]; then
  echo ""
  echo "CHROM unequal"
  exit 1
fi

diff truePOS tmpPOS3
if [ $? -ne 0 ]; then
  echo ""
  echo "POS unequal"
  exit 1
fi


rm tmpCHROM tmpPOS trueCHROM truePOS
