.. _sec-synopsis:

========
Synopsis
========

dEploid [ -vcf *file* ] [ -plaf *file* ] [ -noPanel ] ... [ -exclude *file* ] [ -vcfOut ] [ -o *string* ] \
    \

dEploid [ -vcf *file* ] [ -plaf *file* ] [ -panel *file* ] ... [ -exclude *file* ] [ -vcfOut ] [ -o *string* ] \
    \

.. dEploid [ -ref *file* ] [ -alt *file* ] [ -plaf *file* ]  [ -noPanel ] \
..     \

.. dEploid [ -ref *file* ] [ -alt *file* ] [ -plaf *file* ] [ -panel *file* ] \
..     \


Example:
--------

.. code-block:: bash

    $ ./dEploid -vcf data/exampleData/PG0390-C.eg.vcf.gz \
    -plaf data/exampleData/labStrains.eg.PLAF.txt \
    -exclude data/testData/labStrains.test.exclude.txt \
    -noPanel \
    -o test_run \
    -vcfOut -z
