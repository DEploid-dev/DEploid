PFDECONVVERSION = $(shell git show HEAD | head -1 | sed -e "s/commit //g" | cat)
EXTRA_DIST = bootstrap
COMPILEDATE = $(shell date -u | sed -e "s/ /-/g")
distdir = $(PACKAGE)-$(VERSION)

bin_PROGRAMS = pfDeconv pfDeconv_dbg

TESTS = unit_tests
check_PROGRAMS = unit_tests pfDeconv_dbg pfDeconv_prof
PROG = PFDECONV

common_flags = -std=c++0x -Isrc/ -Isrc/codeCogs/ -Isrc/random/ -DPFDECONVVERSION=\"${PFDECONVVERSION}\" -DCOMPILEDATE=\"${COMPILEDATE}\"

common_src = src/random/fastfunc.cpp \
			 src/random/random_generator.cpp \
			 src/random/mersenne_twister.cpp \
			 src/mcmc.cpp \
			 src/panel.cpp \
			 src/pfDeconvIO.cpp \
			 src/utility.cpp \
			 src/updateHap.cpp \
			 src/atMarker.cpp \
			 src/export/pfDeconvIOExport.cpp \
			 src/export/pfDeconvIOExportPosteriorProb.cpp

debug_src = src/debug/mcmcDebug.cpp

pfDeconv_SOURCES = src/pfDeconv.cpp $(common_src)
pfDeconv_dbg_SOURCES =  $(debug_src) $(pfDeconv_SOURCES)
pfDeconv_prof_SOURCES = $(pfDeconv_SOURCES)

pfDeconv_CXXFLAGS = $(common_flags) -DNDEBUG -O3
pfDeconv_dbg_CXXFLAGS = -g $(common_flags) -O3
pfDeconv_prof_CXXFLAGS = $(common_flags) -DNDEBUG -fno-omit-frame-pointer -pg -O1


unit_tests_SOURCES = $(common_src) \
					 tests/unittest/test_panel.cpp \
					 tests/unittest/test_pfDeconvIO.cpp \
					 tests/unittest/test_runner.cpp \
					 tests/unittest/test_updatePairHap.cpp \
					 tests/unittest/test_updateSingleHap.cpp \
					 tests/unittest/test_atMarker.cpp \
					 tests/unittest/test_utilities.cpp \
					 tests/unittest/test_mcmc.cpp

unit_tests_CXXFLAGS = $(pfDeconv_CXXFLAGS) -DUNITTEST -Wno-write-strings --coverage
unit_tests_LDADD          = -lcppunit -ldl

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf src/*.gc* tests/unittest/*.gc* src/*.gcov
