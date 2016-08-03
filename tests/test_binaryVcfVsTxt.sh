#!/bin/bash
./dEploid -vcf tests/testData/PG0389-C100.vcf  -plaf tests/testData/plafForTesting.csv -panel tests/testData/panelForTesting.csv -o vcf -seed 1
./dEploid -ref tests/testData/PG0389.C_ref100.txt -alt tests/testData/PG0389.C_alt100.txt -plaf tests/testData/plafForTesting.csv -panel tests/testData/panelForTesting.csv -o txt -seed 1

diff txt.prop vcf.prop
if [ $? -ne 0 ]; then
  echo ""
  echo "Proportion unequal"
  exit 1
fi

diff txt.llk vcf.llk
if [ $? -ne 0 ]; then
  echo ""
  echo "Likelihood unequal"
  exit 1
fi

diff txt.hap vcf.hap
if [ $? -ne 0 ]; then
  echo ""
  echo "Haplotypes unequal"
  exit 1
fi
