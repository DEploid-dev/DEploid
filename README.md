PfDeconv
===========

[![Circle CI](https://circleci.com/gh/mcveanlab/PfDeconv.svg?style=svg)](https://circleci.com/gh/mcveanlab/PfDeconv)
[![Build Status](https://travis-ci.org/mcveanlab/PfDeconv.svg?branch=master)](https://travis-ci.org/mcveanlab/PfDeconv)
[![Coverage Status](https://coveralls.io/repos/github/mcveanlab/PfDeconv/badge.svg?branch=coverage)](https://coveralls.io/github/mcveanlab/PfDeconv?branch=coverage)

_PfDeconv_ is developed as part of the [_Pf3k_](https://www.malariagen.net/projects/parasite/pf3k) project. The _Pf3k_ project is a global collaboration using the latest sequencing technologies to provide a high-resolution view of natural variation in the malaria parasite Plasmodium falciparum. Parasite DNA were extracted from patient blood sample, which often contains more than one parasite strain, with unknown proportions. _PfDeconv_ is used for deconvoluting mixed haplotypes, and reporting the mixture proportions from each sample.

##NEWS
(to come)

##INSTALL
```bash
./bootstrap
make
```

##LICENCE
You can freely use all code in this project under the conditions of the GNU GPL Version 3 or later.

##HOW IT WORKS

Program parameters and options:

Options              | Useage |
:-------------------:| ------------------------------- |
-h or -help          |  Help. List the following content.
            -ref STR |  File path of reference allele count.
            -alt STR |  File path of alternative allele count.
           -plaf STR |  File path of population level allele frequencies.
          -panel STR |  File path of the reference panel.
              -o STR |  Specify the file name prefix of the output.
              -p INT |  Out put precision (default value 8).
              -k INT |  Number of strain (default value 5).
           -seed INT |  Random seed.
        -nSample INT |  Number of MCMC samples.
           -rate INT |  MCMC sample rate.
            -noPanel |  Use population level allele frequency as prior.

##Examples:
```bash
./pfDeconv -ref labStrains/PG0390_first100ref.txt -alt labStrains/PG0390_first100alt.txt -plaf labStrains/labStrains_first100_PLAF.txt -panel labStrains/lab_first100_Panel.txt -o tmp1
./pfDeconv -ref labStrains/PG0390_first100ref.txt -alt labStrains/PG0390_first100alt.txt -plaf labStrains/labStrains_first100_PLAF.txt -panel labStrains/lab_first100_Panel.txt -nSample 100 -rate 3
```
