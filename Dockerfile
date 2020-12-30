FROM ubuntu:18.04
MAINTAINER Joe Zhu <sha.joe.zhu@gmail.com>
RUN apt-get update -qq \
    && apt-get install -qq git build-essential autoconf autoconf-archive libcppunit-dev zlib1g-dev bash-completion \
    && apt-cache policy zlib*
RUN git clone --recursive https://github.com/DEploid-dev/DEploid.git
RUN cd /DEploid/ \
    && ./bootstrap \
    && make install
ENTRYPOINT ["dEploid"]
