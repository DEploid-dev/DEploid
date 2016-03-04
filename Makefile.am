PFDECONVVERSION = $(shell git show HEAD | head -1 | sed -e "s/commit //g" | cat)
EXTRA_DIST = bootstrap
COMPILEDATE = $(shell date -u | sed -e "s/ /-/g")
distdir = $(PACKAGE)-$(VERSION)

bin_PROGRAMS = pfDeconv pfDeconv_dbg

TESTS = unit_tests
check_PROGRAMS = unit_tests pfDeconv_dbg pfDeconv_prof
PROG = PFDECONV

common_flags = -std=c++0x -Isrc/ -Isrc/codeCogs/ -DPFDECONVVERSION=\"${PFDECONVVERSION}\" -DCOMPILEDATE=\"${COMPILEDATE}\"

common_src = src/fastfunc.cpp \
			 src/random_generator.cpp \
			 src/mersenne_twister.cpp \
			 src/mcmc.cpp \
			 src/panel.cpp \
			 src/pfDeconvIO.cpp \
			 src/pfDeconvIOExport.cpp \
			 src/utility.cpp \
			 src/updateHap.cpp \
			 src/atMarker.cpp

pfDeconv_SOURCES = src/pfDeconv.cpp $(common_src)
pfDeconv_dbg_SOURCES = $(pfDeconv_SOURCES)
pfDeconv_prof_SOURCES = $(pfDeconv_SOURCES)

pfDeconv_CXXFLAGS = $(common_flags) -DNDEBUG -O3
pfDeconv_dbg_CXXFLAGS = -g $(common_flags) -O3
pfDeconv_prof_CXXFLAGS = $(common_flags) -DNDEBUG -fno-omit-frame-pointer -pg -O1


unit_tests_SOURCES = $(common_src) \
					 tests/unittest/test_panel.cpp \
					 tests/unittest/test_pfDeconvIO.cpp \
					 tests/unittest/test_runner.cpp \
					 tests/unittest/test_updatePairHap.cpp

unit_tests_CXXFLAGS = $(pfDeconv_CXXFLAGS) -Wno-write-strings
unit_tests_LDADD          = -lcppunit -ldl

