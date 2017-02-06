#!/bin/bash
./dEploid -ref data/testData/PG0390-C.test.ref \
-alt data/testData/PG0390-C.test.alt \
-plaf data/testData/labStrains.test.PLAF.txt \
-panel data/testData/labStrains.test.panel.txt -o tmp -exportPostProb -k 2
cut -f 1 tmp.hap | tail -n+2 > tmpCHROM
cut -f 2 tmp.hap | tail -n+2 > tmpPOS


cut -f 1 data/testData/PG0390-C.test.ref  | tail -n+2 > trueCHROM
cut -f 2 data/testData/PG0390-C.test.ref | tail -n+2 > truePOS


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

cut -f 1 tmp.single0 | tail -n+2 > tmpCHROM0
cut -f 2 tmp.single0 | tail -n+2 > tmpPOS0

diff trueCHROM tmpCHROM0
if [ $? -ne 0 ]; then
  echo ""
  echo "CHROM unequal"
  exit 1
fi

diff truePOS tmpPOS0
if [ $? -ne 0 ]; then
  echo ""
  echo "POS unequal"
  exit 1
fi


#cut -f 1 tmp.pair | tail -n+2 > tmpCHROM3
#cut -f 2 tmp.pair | tail -n+2 > tmpPOS3

#diff trueCHROM tmpCHROM3
#if [ $? -ne 0 ]; then
  #echo ""
  #echo "CHROM unequal"
  #exit 1
#fi

#diff truePOS tmpPOS3
#if [ $? -ne 0 ]; then
  #echo ""
  #echo "POS unequal"
  #exit 1
#fi


rm tmpCHROM tmpPOS trueCHROM truePOS
