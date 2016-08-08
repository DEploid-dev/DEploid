=======
DEploid
=======

[![Circle CI](https://circleci.com/gh/mcveanlab/DEploid.svg?style=svg)](https://circleci.com/gh/mcveanlab/DEploid)
[![Build Status](https://travis-ci.org/mcveanlab/DEploid.svg?branch=master)](https://travis-ci.org/mcveanlab/DEploid)
[![Coverage Status](https://coveralls.io/repos/github/mcveanlab/DEploid/badge.svg)](https://coveralls.io/github/mcveanlab/DEploid)


************
Installation
************

    ./bootstrap
    make


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


*******
Licence
*******
You can freely use all code in this project under the conditions of the GNU GPL Version 3 or later.





