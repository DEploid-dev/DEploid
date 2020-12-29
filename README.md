
<img src="docs/_static/deploidTrans.png" width="180">


[![License (GPL version 3)](https://img.shields.io/badge/license-GPL%20version%203-brightgreen.svg)](http://opensource.org/licenses/GPL-3.0)
[![Build Status](https://travis-ci.org/DEploid-dev/DEploid.svg?branch=master)](https://travis-ci.org/DEploid-dev/DEploid)
[![CircleCI](https://circleci.com/gh/DEploid-dev/DEploid.svg?style=shield)](https://circleci.com/gh/DEploid-dev/DEploid)
[![Coverage Status](https://coveralls.io/repos/github/DEploid-dev/DEploid/badge.svg)](https://coveralls.io/github/DEploid-dev/DEploid)
[![Documentation Status](http://readthedocs.org/projects/deploid/badge/?version=latest)](http://deploid.readthedocs.io/en/latest/)
[![Docker Status](https://img.shields.io/docker/build/shajoezhu/deploid.svg)](https://hub.docker.com/r/shajoezhu/deploid/)

`dEploid` is designed for deconvoluting mixed genomes with unknown proportions. Traditional ‘phasing’ programs are limited to diploid organisms. Our method modifies Li and Stephen’s algorithm with Markov chain Monte Carlo (MCMC) approaches, and builds a generic framework that allows haloptype searches in a multiple infection setting.

Please see the [documentation](http://deploid.readthedocs.io/en/latest/) for further details.

Installation
------------

You can also install `dEploid` directly from the git repository. Here, you will need `autoconf`, check whether this is already installed by running:
```bash
$ which autoconf
```

On Debian/Ubuntu based systems:
```bash
$ apt-get install build-essential autoconf autoconf-archive libcppunit-dev zlib1g-dev
```

On Mac OS:
```bash
$ port install automake autoconf autoconf-archive cppunit
```

Afterwards you can clone the code from the github repository,
```bash
$ git clone git@github.com:mcveanlab/DEploid.git
$ cd DEploid
$ git submodule update --init --recursive --remote
```

and build the binary using
```bash
$ ./bootstrap
$ make
```

Usage
-----

Please see the [documentation](http://deploid.readthedocs.io/en/latest/) for further details.


Docker image
------------

```bash
docker pull shajoezhu/deploid
docker run -v ${PWD}:/tmp/ -w /tmp/  shajoezhu/deploid ...
```

Licence
-------

You can freely use all code in this project under the conditions of the GNU GPL Version 3 or later.


Citation
--------

If you use `dEploid` with the flag `-ibd`, please cite the following paper:

Zhu, J. S., J. A. Hendry, J. Almagro-Garcia, R. D. Pearson, R. Amato, A. Miles, D. J. Weiss, T. C. D. Lucas, M. Nguyen, P. W. Gething, D. Kwiatkowski, G. McVean, and for the Pf3k Project. (2018) The origins and relatedness structure of mixed infections vary with local prevalence of *P. falciparum* malaria. *eLife*, 40845, doi: https://doi.org/10.7554/eLife.40845.


If you use `dEploid` in your work, please cite the program:

Zhu, J. S., J. A. Garcia, G. McVean. (2018) Deconvolution of multiple infections in *Plasmodium falciparum* from high throughput sequencing data. *Bioinformatics* 34(1), 9-15. doi: https://doi.org/10.1093/bioinformatics/btx530.
