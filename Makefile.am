DEPLOIDVERSION = $(shell git show HEAD | head -1 | sed -e "s/commit //g" | cat)
EXTRA_DIST = bootstrap
COMPILEDATE = $(shell date -u | sed -e "s/ /-/g")
distdir = $(PACKAGE)-$(VERSION)

bin_PROGRAMS = dEploid dEploid_dbg

man1_MANS = docs/_build/man/dEploid.1

TESTS = unit_tests
check_PROGRAMS = unit_tests dEploid_dbg dEploid_prof
PROG = DEPLOID

common_flags = -std=c++0x -Isrc/ -Isrc/codeCogs/ -Isrc/random/ -DDEPLOIDVERSION=\"${DEPLOIDVERSION}\" -DCOMPILEDATE=\"${COMPILEDATE}\"

common_src = src/random/fastfunc.cpp \
			 src/random/random_generator.cpp \
			 src/random/mersenne_twister.cpp \
			 src/mcmc.cpp \
			 src/panel.cpp \
			 src/dEploidIO.cpp \
			 src/utility.cpp \
			 src/updateHap.cpp \
			 src/atMarker.cpp \
			 src/export/dEploidIOExport.cpp \
			 src/export/dEploidIOExportPosteriorProb.cpp \
			 src/vcfReader.cpp

debug_src = src/debug/mcmcDebug.cpp src/debug/vcfReaderDebug.cpp

dEploid_SOURCES = src/dEploid.cpp $(common_src)
dEploid_dbg_SOURCES =  $(debug_src) $(dEploid_SOURCES)
dEploid_prof_SOURCES = $(dEploid_SOURCES)

dEploid_CXXFLAGS = $(common_flags) -DNDEBUG -O3
dEploid_dbg_CXXFLAGS = -g $(common_flags) -O3
dEploid_prof_CXXFLAGS = $(common_flags) -DNDEBUG -fno-omit-frame-pointer -pg -O1


unit_tests_SOURCES = $(common_src) \
					 tests/unittest/test_panel.cpp \
					 tests/unittest/test_dEploidIO.cpp \
					 tests/unittest/test_runner.cpp \
					 tests/unittest/test_updatePairHap.cpp \
					 tests/unittest/test_updateSingleHap.cpp \
					 tests/unittest/test_atMarker.cpp \
					 tests/unittest/test_utilities.cpp \
					 tests/unittest/test_mcmc.cpp \
					 tests/unittest/test_atMarker.cpp \
					 tests/unittest/test_vcfReader.cpp

unit_tests_CXXFLAGS = $(common_flags) -DNDEBUG -DUNITTEST -Wno-write-strings --coverage
unit_tests_LDADD          = -lcppunit -ldl

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf src/*.gc* tests/unittest/*.gc* src/*.gcov
