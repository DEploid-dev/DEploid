#!/bin/bash
sameFlags="-exclude data/testData/labStrains.test.exclude.txt -plaf data/testData/labStrains.test.PLAF.txt -panel data/testData/labStrains.test.panel.txt -seed 1 -vcfOut -vcf data/testData/PG0390-C.test.vcf.gz"
./dEploid ${sameFlags} -o vcf1 -z || exit 1
./dEploid ${sameFlags} -o vcf2 -z || exit 1

jobBrief="classic"

gzip -d vcf1.${jobBrief}.vcf.gz
gzip -d vcf2.${jobBrief}.vcf.gz

diff vcf2.${jobBrief}.prop vcf1.${jobBrief}.prop
if [ $? -ne 0 ]; then
  echo ""
  echo "Proportion unequal"
  exit 1
fi

diff vcf2.${jobBrief}.llk vcf1.${jobBrief}.llk
if [ $? -ne 0 ]; then
  echo ""
  echo "Likelihood unequal"
  exit 1
fi

diff vcf2.${jobBrief}.hap vcf1.${jobBrief}.hap
if [ $? -ne 0 ]; then
  echo ""
  echo "Haplotypes unequal"
  exit 1
fi

grep "##" data/testData/PG0390-C.test.vcf > originalHeader
head -165 vcf1.${jobBrief}.vcf > newHeader
diff originalHeader newHeader
if [ $? -ne 0 ]; then
  echo ""
  echo "Vcf header unequal"
  exit 1
fi

# Because of exclude, do not compare vcf columns 1 to 9

tail -n +2 vcf2.${jobBrief}.hap | cut -f 1-2 > txtHap1to2
grep -v "#" vcf2.${jobBrief}.vcf | cut -f 1-2 > txtVcf1to2
diff txtHap1to2 txtVcf1to2
if [ $? -ne 0 ]; then
  echo ""
  echo "Text Vcf columns unequal"
  exit 1
fi

tail -n +2 vcf2.${jobBrief}.hap | cut -f 3- > txtHap
grep -v "#" vcf1.${jobBrief}.vcf | cut -f 10- > newVcfHap
diff txtHap newVcfHap
if [ $? -ne 0 ]; then
  echo ""
  echo "Vcf Hap unequal"
  exit 1
fi

grep -v "#" vcf2.${jobBrief}.vcf | cut -f 10- > newTxtHap
diff txtHap newTxtHap
if [ $? -ne 0 ]; then
  echo ""
  echo "Text Vcf columns unequal"
  exit 1
fi

rm vcf* txt* newHeader newTxtHap newVcfHap originalHeader
