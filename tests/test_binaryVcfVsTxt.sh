#!/bin/bash
#./dEploid -vcf tests/testData/PG0389-C100.vcf  -plaf tests/testData/plafForTesting.csv -panel tests/testData/panelForTesting.csv -o vcf -seed 1 -vcfOut
#./dEploid -ref tests/testData/PG0389.C_ref100.txt -alt tests/testData/PG0389.C_alt100.txt -plaf tests/testData/plafForTesting.csv -panel tests/testData/panelForTesting.csv -o txt -seed 1  -vcfOut

./dEploid -vcf data/testData/PG0390-C.test.vcf -exclude data/testData/labStrains.test.exclude.txt -plaf data/testData/labStrains.test.PLAF.txt -panel data/testData/labStrains.test.panel.txt -o vcf -seed 1 -vcfOut || exit 1
./dEploid -ref data/testData/PG0390-C.test.ref -alt data/testData/PG0390-C.test.alt -exclude data/testData/labStrains.test.exclude.txt -plaf data/testData/labStrains.test.PLAF.txt -panel data/testData/labStrains.test.panel.txt -o txt -seed 1  -vcfOut || exit 1

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

grep "##" data/testData/PG0390-C.test.vcf > originalHeader
head -165 vcf.vcf > newHeader
diff originalHeader newHeader
if [ $? -ne 0 ]; then
  echo ""
  echo "Vcf header unequal"
  exit 1
fi

grep -v "#" data/testData/PG0390-C.test.vcf | cut -f 1-8 > originalVcf1to8
grep -v "#" vcf.vcf | cut -f 1-8 > newVcf1to8
diff originalVcf1to8 newVcf1to8
if [ $? -ne 0 ]; then
  echo ""
  echo "Vcf columns unequal"
  exit 1
fi

tail -n +2 txt.hap | cut -f 1-2 > txtHap1to2
grep -v "#" txt.vcf | cut -f 1-2 > txtVcf1to2
diff txtHap1to2 txtVcf1to2
if [ $? -ne 0 ]; then
  echo ""
  echo "Text Vcf columns unequal"
  exit 1
fi

tail -n +2 txt.hap | cut -f 3- > txtHap
grep -v "#" vcf.vcf | cut -f 10- > newVcfHap
diff txtHap newVcfHap
if [ $? -ne 0 ]; then
  echo ""
  echo "Vcf Hap unequal"
  exit 1
fi

grep -v "#" txt.vcf | cut -f 10- > newTxtHap
diff txtHap newTxtHap
if [ $? -ne 0 ]; then
  echo ""
  echo "Text Vcf columns unequal"
  exit 1
fi

rm vcf* txt*
