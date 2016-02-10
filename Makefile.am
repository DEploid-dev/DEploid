#distdir = $(PACKAGE)-$(VERSION)

bin_PROGRAMS = pfDeconv

TESTS = unit_tests
check_PROGRAMS = unit_tests
PROG = PFDECONV

common_flags = -DNDEBUG -O3 -g -std=c++0x -Isrc/

#COMPILEDATE = $(shell date -u | sed -e "s/ /-/g")

pfDeconv_SOURCES = src/pfDeconv.cpp \
				 src/fastfunc.cpp \
				 src/random_generator.cpp \
				 src/mersenne_twister.cpp


pfDeconv_CXXFLAGS = $(common_flags) -DNDEBUG

unit_tests_SOURCES = tests/unittest/test_runner.cpp
unit_tests_CXXFLAGS = $(pfDeconv_CXXFLAGS)
unit_tests_LDADD          = -lcppunit -ldl

