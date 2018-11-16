#!/bin/bash
sameFlags="-exclude data/testData/labStrains.test.exclude.txt -plaf data/testData/labStrains.test.PLAF.txt -panel data/testData/labStrains.test.panel.txt -seed 1 -vcfOut "
./dEploid ${sameFlags} -vcf data/testData/PG0390-C.test.vcf.gz -o vcf -z || exit 1
./dEploid ${sameFlags} -ref data/testData/PG0390-C.test.ref -alt data/testData/PG0390-C.test.alt -o txt || exit 1

jobBrief="classic"

gzip -d vcf.${jobBrief}.vcf.gz


diff txt.${jobBrief}.prop vcf.${jobBrief}.prop
if [ $? -ne 0 ]; then
  echo ""
  echo "Proportion unequal"
  exit 1
fi

diff txt.${jobBrief}.llk vcf.${jobBrief}.llk
if [ $? -ne 0 ]; then
  echo ""
  echo "Likelihood unequal"
  exit 1
fi

diff txt.${jobBrief}.hap vcf.${jobBrief}.hap
if [ $? -ne 0 ]; then
  echo ""
  echo "Haplotypes unequal"
  exit 1
fi

grep "##" data/testData/PG0390-C.test.vcf > originalHeader
head -165 vcf.${jobBrief}.vcf > newHeader
diff originalHeader newHeader
if [ $? -ne 0 ]; then
  echo ""
  echo "Vcf header unequal"
  exit 1
fi

# Because of exclude, do not compare vcf columns 1 to 9

tail -n +2 txt.${jobBrief}.hap | cut -f 1-2 > txtHap1to2
grep -v "#" txt.${jobBrief}.vcf | cut -f 1-2 > txtVcf1to2
diff txtHap1to2 txtVcf1to2
if [ $? -ne 0 ]; then
  echo ""
  echo "Text Vcf columns unequal"
  exit 1
fi

tail -n +2 txt.${jobBrief}.hap | cut -f 3- > txtHap
grep -v "#" vcf.${jobBrief}.vcf | cut -f 10- > newVcfHap
diff txtHap newVcfHap
if [ $? -ne 0 ]; then
  echo ""
  echo "Vcf Hap unequal"
  exit 1
fi

grep -v "#" txt.${jobBrief}.vcf | cut -f 10- > newTxtHap
diff txtHap newTxtHap
if [ $? -ne 0 ]; then
  echo ""
  echo "Text Vcf columns unequal"
  exit 1
fi

rm vcf* txt*
