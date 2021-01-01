Synopsis
========

dEploid [ -vcf *file* ] [ -sample *string* ] [ -plafFromVcf ] [ -panel *file* ] ... [ -o *string* ]


Example:
--------

```bash
    $ ./dEploid -vcf data/testData/PG0390-C.test.vcf.gz \
    -sample PG0390-C -plafFromVcf \
    -exclude data/testData/labStrains.test.exclude.txt.gz \
    -panel data/testData/labStrains.test.panel.txt.gz \
    -o test_run \
    -vcfOut -z \
    -best
```

For help and more examples, type `dEploid -help`.
