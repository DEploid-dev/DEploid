This directory contains test data which are used for unit testing and testing for binary

Data files include crappy vcf files at `tests/testData/crappyVcf`, which are used for testing bad vcf exceptions.

Files for testing the txtReader exclude site function is working
```
data/testData/txtReaderForTestingAfterExclude.txt
data/testData/txtReaderForTestingToBeExclude.txt
data/testData/txtReaderForTesting.txt
```
exclude markers in `data/testData/txtReaderForTesting.txt` from list `data/testData/txtReaderForTestingToBeExclude.txt`, should get `data/testData/txtReaderForTestingAfterExclude.txt`

All the rest of the test files focus on pf3k sample PG0390-C

