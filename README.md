
<img src="docs/_static/deploid.png" width="180">
[![Circle CI](https://circleci.com/gh/mcveanlab/DEploid.svg?style=svg)](https://circleci.com/gh/mcveanlab/DEploid)
[![Build Status](https://travis-ci.org/mcveanlab/DEploid.svg?branch=master)](https://travis-ci.org/mcveanlab/DEploid)
[![Coverage Status](https://coveralls.io/repos/github/mcveanlab/DEploid/badge.svg)](https://coveralls.io/github/mcveanlab/DEploid)

`dEploid` is designed for deconvoluting mixed genomes with unknown proportions. Traditional ‘phasing’ programs are limited to diploid organisms. Our method modifies Li and Stephen’s algorithm with Markov chain Monte Carlo (MCMC) approaches, and builds a generic framework that allows haloptype searches in a multiple infection setting.

Please see the [documentation](http://deploid.readthedocs.io/en/latest/) for further details.

Installation
------------

You can also install `dEploid` directly from the git repository. Here, you need to install `autoconf` first:

On Debian/Ubuntu based systems:
```bash
apt-get install build-essential autoconf autoconf-archive libcppunit-dev
```

On Mac OS:
```bash
port install automake autoconf autoconf-archive cppunit
```

Afterwards you can build the binary using
```bash
./bootstrap
make
```

Usage
-----

Please see the [documentation](http://deploid.readthedocs.io/en/latest/) for further details.


Licence
-------

You can freely use all code in this project under the conditions of the GNU GPL Version 3 or later.


Citation
--------

If you use `dEploid` in your work, please cite the program:

PLACEHOLDER FOR APP NOTE



