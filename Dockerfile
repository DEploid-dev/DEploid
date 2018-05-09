MAINTAINER Joe Zhu <joe.zhu@bdi.ox.ac.uk>

RUN git clone --recursive https://github.com/mcveanlab/DEploid.git
WORKDIR /DEploid
RUN ./bootstrap
RUN make
CMD ./dEploid

#  - if [ $TRAVIS_OS_NAME == linux ]; then sudo apt-get update -qq; sudo apt-get install -qq libcppunit-dev valgrind r-base-core lcov python-pip doxygen graphviz; pip install --user cpp-coveralls; fi
#  - if [ $TRAVIS_OS_NAME == linux ]; then apt-cache policy zlib*; fi

