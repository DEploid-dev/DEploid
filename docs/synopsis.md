Synopsis
========

dEploid [ -vcf *file* ] [ -plaf *file* ] [ -noPanel ] ... [ -exclude *file* ] [ -vcfOut ] [ -o *string* ]


dEploid [ -vcf *file* ] [ -plaf *file* ] [ -panel *file* ] ... [ -exclude *file* ] [ -vcfOut ] [ -o *string* ]


Example:
--------

```bash
$ ./dEploid -vcf data/exampleData/PG0390-C.eg.vcf.gz \
-plaf data/exampleData/labStrains.eg.PLAF.txt \
-exclude data/testData/labStrains.test.exclude.txt \
-noPanel \
-o test_run \
-vcfOut -z
```
