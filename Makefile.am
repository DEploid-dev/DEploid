#PFDECONVVERSION = $(shell git show HEAD | head -1 | sed -e "s/commit //g" | cat)
EXTRA_DIST = bootstrap
#COMPILEDATE = $(shell date -u | sed -e "s/ /-/g")
#distdir = $(PACKAGE)-$(VERSION)

bin_PROGRAMS = pfDeconv

TESTS = unit_tests
check_PROGRAMS = unit_tests pfDeconv_dbg
PROG = PFDECONV

#common_flags = -O3 -std=c++0x -Isrc/ -DPFDECONVVERSION=\"${PFDECONVVERSION}\" -DCOMPILEDATE=\"${COMPILEDATE}\"

common_flags = -O3 -std=c++0x -Isrc/

pfDeconv_SOURCES = src/pfDeconv.cpp \
				 src/fastfunc.cpp \
				 src/random_generator.cpp \
				 src/mersenne_twister.cpp \
				 src/mcmc.cpp

pfDeconv_dbg_SOURCES = $(pfDeconv_SOURCES)

pfDeconv_CXXFLAGS = $(common_flags) -DNDEBUG
pfDeconv_dbg_CXXFLAGS = -g $(common_flags)

unit_tests_SOURCES = tests/unittest/test_runner.cpp
unit_tests_CXXFLAGS = $(pfDeconv_CXXFLAGS)
unit_tests_LDADD          = -lcppunit -ldl

