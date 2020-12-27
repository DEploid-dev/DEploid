/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2018 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of dEploid.
 *
 * dEploid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "src/dEploidIO.hpp"

class TestIO : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE(TestIO);
    CPPUNIT_TEST(testMainConstructor);
    CPPUNIT_TEST(testInitialization);
    CPPUNIT_TEST(testPrintHelp);
    CPPUNIT_TEST(testPrintVersion);
    CPPUNIT_TEST(testNotEnoughArg);
    CPPUNIT_TEST(testOutOfRange);
    CPPUNIT_TEST(testWrongType);
    CPPUNIT_TEST(testFileNameMissing);
    CPPUNIT_TEST(testInvalidInputFile);
    CPPUNIT_TEST(testUnknowArg);
    CPPUNIT_TEST(testFlagsConflict);
    CPPUNIT_TEST(testExtractRefAltPlaf);
    CPPUNIT_TEST(testLociNumberUnequal);
    CPPUNIT_TEST(testForbidMoves);
    CPPUNIT_TEST(testInitialProp);
    CPPUNIT_TEST(testInvalidVcf);
    CPPUNIT_TEST(testVcfHeader);
    CPPUNIT_TEST(testVcfNoAD);
    CPPUNIT_TEST(testInvalidVcfGz);
    CPPUNIT_TEST(testVcfGzHeader);
    CPPUNIT_TEST(testVcfGzNoAD);
    CPPUNIT_TEST(testVcfOutUnSpecified);
    CPPUNIT_TEST(testVcfGzNoVQSLOD);
    CPPUNIT_TEST(testInvalidK);
    CPPUNIT_TEST(testChromPainting);
    CPPUNIT_TEST(testComputeObsWsaf);
    CPPUNIT_TEST_SUITE_END();

 private:
    double epsilon2;
    double epsilon3;

 public:
    void setUp() {
        this->epsilon3 = 0.000000000001;
        this->epsilon2 = 0.0000001;
    }

    void tearDown() {
    }

    void testInitialization() {
        DEploidIO* dEploidIOptr = new DEploidIO();
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->randomSeedWasGiven(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->doExportRecombProb(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->compressVcf(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->initialPropWasGiven(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->excludeSites(), false);
        CPPUNIT_ASSERT(dEploidIOptr->excludedMarkers == NULL);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->randomSeed(), (size_t)0);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->help(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->usePanel(), true);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->precision_, (size_t)8);
        CPPUNIT_ASSERT(dEploidIOptr->prefix_ == "pf3k-dEploid");
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->kStrain_, (size_t)4);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->nMcmcSample_, (size_t)800);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->doUpdateProp(), true);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->doUpdatePair(), true);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->doUpdateSingle(), true);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->doExportPostProb(), false);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(dEploidIOptr->mcmcBurn_, 0.5, epsilon3);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->mcmcMachineryRate_, (size_t)5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(dEploidIOptr->missCopyProb_,
                                     0.01, epsilon3);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->useConstRecomb(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->forbidCopyFromSame(), false);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(dEploidIOptr->constRecombProb_,
                                     1.0, epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(dEploidIOptr->averageCentimorganDistance_,
                                     15000.0, epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(dEploidIOptr->parameterG(),
                                     20.0, epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(dEploidIOptr->parameterSigma(),
                                     5.0, epsilon3);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->doExportVcf(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIOptr->doLsPainting(), false);
        delete dEploidIOptr;
    }


    void testMainConstructor() {
        char *argv1[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1",
                         "-exclude",
                         "data/testData/labStrains.test.exclude.txt",
                         "-vcfOut", "-z"};
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(11, argv1));
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(13, argv1));
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(14, argv1));
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(15, argv1));

        char *argv2[] = { "./dEploid",
                         "-vcf", "data/testData/PG0390-C.test.vcf",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1",
                         "-initialHap",
                         "data/testData/PG0390-C.test.nopanel.hap",
                         "-exclude",
                         "data/testData/labStrains.test.exclude.txt"};
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(9, argv2));
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(11, argv2));
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(13, argv2));
    }


    void testConstructorCopier() {
        char *argv1[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1",
                         "-exclude",
                         "data/testData/labStrains.test.exclude.txt",
                         "-vcfOut", "-z"};
        DEploidIO tmp(11, argv1);
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(tmp));
    }

    void testInitialProp() {
        char *argv1[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-initialP", "0.1", "0.2", "0.3", "0.4", "-k", "4" };
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(16, argv1));
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(14, argv1));
        DEploidIO tmp(14, argv1);
        CPPUNIT_ASSERT_EQUAL(tmp.kStrain(), (size_t)4);

        char *argv2[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-initialP", "0.1", "0.2", "0.3", "0.3", "0.2"};
        CPPUNIT_ASSERT_THROW(DEploidIO(15, argv2), SumOfPropNotOne);

        char *argv3[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-initialP", "0.1", "0.2", "0.3", "0.3", "0.02"};
        CPPUNIT_ASSERT_THROW(DEploidIO(15, argv3), SumOfPropNotOne);

        char *argv4[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-initialP", "-o", "tmp"};
        CPPUNIT_ASSERT_THROW(DEploidIO(12, argv4), NotEnoughArg);

        char *argv5[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-initialP"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv5), NotEnoughArg);

        char *argv6[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-initialP", "0.1", "0.2", "0.3", "0.4", "-k", "5" };
        CPPUNIT_ASSERT_THROW(DEploidIO(16, argv6), NumOfPropNotMatchNumStrain);

        char *argv7[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                          "-k", "5", "-initialP", "0.1", "0.2", "0.3", "0.4" };
        CPPUNIT_ASSERT_THROW(DEploidIO(16, argv7), NumOfPropNotMatchNumStrain);
    }


    void testPrintHelp() {
        char *argv[] = { "./dEploid" };
        std::ostream *output = &std::cout;
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(1, argv));
        DEploidIO dEploidIO1(1, argv);
        CPPUNIT_ASSERT_EQUAL((size_t)0, dEploidIO1.argv_.size());
        CPPUNIT_ASSERT_EQUAL(true, dEploidIO1.help());
        CPPUNIT_ASSERT_NO_THROW(dEploidIO1.printHelp(*output));

        char *argv1[] = { "./dEploid", "-h" };
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(2, argv1));
        DEploidIO dEploidIO2(2, argv1);
        CPPUNIT_ASSERT_EQUAL((size_t)1, dEploidIO2.argv_.size());
        CPPUNIT_ASSERT_EQUAL(true, dEploidIO2.help());
        CPPUNIT_ASSERT_NO_THROW(dEploidIO2.printHelp(*output));
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(string("-help")));

        char *argv2[] = { "./dEploid", "-help" };
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(2, argv2));
        DEploidIO dEploidIO3(2, argv2);
        CPPUNIT_ASSERT_EQUAL((size_t)1, dEploidIO3.argv_.size());
        CPPUNIT_ASSERT_EQUAL(true, dEploidIO3.help());
        CPPUNIT_ASSERT_NO_THROW(dEploidIO3.printHelp(*output));
    }


    void testPrintVersion() {
        std::ostream *output = &std::cout;
        char *argv1[] = { "./dEploid", "-v" };
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(2, argv1));
        DEploidIO dEploidIO1(2, argv1);
        CPPUNIT_ASSERT_EQUAL(true, dEploidIO1.version());
        CPPUNIT_ASSERT_NO_THROW(dEploidIO1.printVersion(*output));

        char *argv2[] = { "./dEploid", "-version" };
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(2, argv2));
        DEploidIO dEploidIO2(2, argv2);
        CPPUNIT_ASSERT_EQUAL(true, dEploidIO2.version());
        CPPUNIT_ASSERT_NO_THROW(dEploidIO2.printVersion(*output));
    }


    void testOutOfRange() {
        char *argv1[] = {"./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-burn", "1.1" };
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv1), OutOfRange);

        char *argv2[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-miss", "5.1"};
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv2), OutOfRange);

        char *argv3[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-recomb", "2.1" };
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv3), OutOfRange);
    }


    void testNotEnoughArg() {
        char *argv1[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o" };
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv1), NotEnoughArg);

        char *argv2[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-seed"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv2), NotEnoughArg);

        char *argv3[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-p"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv3), NotEnoughArg);

        char *argv4[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-k"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv4), NotEnoughArg);

        char *argv5[] = { "./dEploid",
                         "-ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt" };
        CPPUNIT_ASSERT_THROW(DEploidIO(8, argv5), NotEnoughArg);

        char *argv6[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt" };
        CPPUNIT_ASSERT_THROW(DEploidIO(8, argv6), NotEnoughArg);

        char *argv7[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf",
                         "-panel", "data/testData/labStrains.test.panel.txt" };
        CPPUNIT_ASSERT_THROW(DEploidIO(8, argv7), NotEnoughArg);

        char *argv8[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel" };
        CPPUNIT_ASSERT_THROW(DEploidIO(8, argv8), NotEnoughArg);

        char *argv9[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-exclude" };
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv9), NotEnoughArg);
    }


    void testWrongType() {
        char *argv2[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-seed", "asdf"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv2), NotEnoughArg);
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv2), WrongType);

        char *argv3[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-p", "asdf"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv3), NotEnoughArg);
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv3), WrongType);

        char *argv4[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-k", "asdf"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv4), NotEnoughArg);
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv4), WrongType);

        char *argv5[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-nSample", "asdf"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv5), NotEnoughArg);
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv5), WrongType);

        char *argv6[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-rate", "asdf"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv6), NotEnoughArg);
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv6), WrongType);
    }


    void testFileNameMissing() {
        char *argv1[] = { "./dEploid",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW(DEploidIO(9, argv1), FileNameMissing);
        CPPUNIT_ASSERT_THROW_MESSAGE(
                "\033[1;31mRef count\033[0m file path missing!",
                DEploidIO(9, argv1), FileNameMissing);
        try {
            DEploidIO(9, argv1);
        }
        catch(const exception &e) {
            CPPUNIT_ASSERT_EQUAL(
                string("\033[1;31mRef count\033[0m file path missing!"),
                string(e.what()));
        }

        char *argv2[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW(DEploidIO(9, argv2), FileNameMissing);

        try {
            DEploidIO(9, argv2);
        }
        catch(const exception &e) {
            CPPUNIT_ASSERT_EQUAL(
               string("\033[1;31mAlt count\033[0m file path missing!"),
               string(e.what()));
        }

        // plaf missing
        char *argv3[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW(DEploidIO(9, argv3), FileNameMissing);

        // panel missing
        char *argv4[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW(DEploidIO(9, argv4), FileNameMissing);
    }


    void testInvalidInputFile() {
        char *argv[] = { "./dEploid",
                         "-ref", "PG0390_first100ref.txt",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv), InvalidInputFile);
        try {
            DEploidIO(11, argv);
        }
        catch(const exception &e) {
          CPPUNIT_ASSERT_EQUAL(
          string("Invalid input file: \033[1;31mPG0390_first100ref.txt\033[0m"),
          string(e.what()));
        }
    }


    void testUnknowArg() {
        char *argv[] = { "./dEploid", "-unknow"};
        CPPUNIT_ASSERT_THROW(DEploidIO(2, argv), UnknowArg);
    }


    void testFlagsConflict() {
        // panel conflict with noPanel
        char *argv1[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-noPanel" };
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv1), FlagsConflict);

        // noPanel conflict with panel
        char *argv2[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-noPanel",
                         "-panel", "data/testData/labStrains.test.panel.txt"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv2), FlagsConflict);

        // ref, alt conflict with vcf
        char *argv3[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1",
                         "-vcf", "data/testData/PG0390-C.test.vcf"};
        CPPUNIT_ASSERT_THROW(DEploidIO(13, argv3), FlagsConflict);

        // vcf conflict with ref
        char *argv4[] = { "./dEploid",
                         "-vcf", "data/testData/PG0390-C.test.vcf",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW(DEploidIO(13, argv4), FlagsConflict);

        // vcf conflict with alt
        char *argv5[] = { "./dEploid",
                         "-vcf", "data/testData/PG0390-C.test.vcf",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW(DEploidIO(13, argv5), FlagsConflict);

        // exportPostProb conflict with noPanel
        char *argv6[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-noPanel", "-exportPostProb" };
        CPPUNIT_ASSERT_THROW(DEploidIO(9, argv6), FlagsConflict);

        // noPanel conflict with exportPostProb
        char *argv7[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-exportPostProb", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(9, argv7), FlagsConflict);

        // best conflict with noPanel
        char *argv8[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-noPanel", "-best" };
        CPPUNIT_ASSERT_THROW(DEploidIO(9, argv8), FlagsConflict);

        // noPanel conflict with best
        char *argv9[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-best", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(9, argv9), FlagsConflict);

        // best conflict with ibd
        char *argv10[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-ibd", "-best", "-noPanel" };
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv10), FlagsConflict);

        // ibd conflict with best
        char *argv11[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-best", "-ibd", "-noPanel" };
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv11), FlagsConflict);

        // best conflict with lasso
        char *argv12[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-lasso", "-best"
                         };
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv12), FlagsConflict);

        // lasso conflict with best
        char *argv13[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-best", "-lasso"
                         };
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv13), FlagsConflict);

        // lasso conflict with ibd
        char *argv14[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-ibd", "-lasso"
                         };
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv14), FlagsConflict);

        // ibd conflict with lasso
        char *argv15[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-lasso", "-ibd"
                         };
        CPPUNIT_ASSERT_THROW(DEploidIO(11, argv15), FlagsConflict);

        // plaf conflict with plafFromVcf
        char *argv16[] = { "./dEploid",
                         "-vcf", "data/testData/PG0390-C.test.vcf",
                         "-plafFromVcf",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt"
                         };
        CPPUNIT_ASSERT_THROW(DEploidIO(8, argv16), FlagsConflict);
    }


    void testLociNumberUnequal() {
        char *argv1[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390.C_alt.txt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt" };
        CPPUNIT_ASSERT_THROW(DEploidIO(9, argv1), LociNumberUnequal);

        char *argv2[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains_samples_PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt" };
        CPPUNIT_ASSERT_THROW(DEploidIO(9, argv2), LociNumberUnequal);
    }


    void testExtractRefAltPlaf() {
        char *argv1[] = { "./dEploid",
         "-ref", "data/testData/PG0390-C.test.ref",
         "-alt", "data/testData/PG0390-C.test.alt",
         "-plaf", "data/testData/labStrains.test.PLAF.txt",
         "-panel", "data/testData/labStrains.test.panel.txt",
         "-o", "tmp1" };


        CPPUNIT_ASSERT_NO_THROW(DEploidIO(11, argv1));
        DEploidIO dEploidIO1(11, argv1);
        CPPUNIT_ASSERT_EQUAL((size_t)594, dEploidIO1.nLoci_);
        // Pf3D7_01_v3 93157 0.35824742268032
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.35824742268032,
                                     dEploidIO1.plaf_[0], 0.000000000001);
        // Pf3D7_01_v3 93157 85
        CPPUNIT_ASSERT_EQUAL(85.0, dEploidIO1.refCount_[0]);
        // Pf3D7_01_v3 93157 0
        CPPUNIT_ASSERT_EQUAL(0.0, dEploidIO1.altCount_[0]);

        // Pf3D7_01_v3 95518 0.553038882597102
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.553038882597102,
                                     dEploidIO1.plaf_[4], 0.000000000001);
        // Pf3D7_01_v3 95518 156
        CPPUNIT_ASSERT_EQUAL(156.0, dEploidIO1.refCount_[4]);
        // Pf3D7_01_v3 95518 46
        CPPUNIT_ASSERT_EQUAL(46.0, dEploidIO1.altCount_[4]);

        // Pf3D7_01_v3 112038 0.178471474703944
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.178471474703944,
                                     dEploidIO1.plaf_[20], 0.000000000001);
        // Pf3D7_01_v3 112038 157
        CPPUNIT_ASSERT_EQUAL(157.0, dEploidIO1.refCount_[20]);
        // Pf3D7_01_v3 112038 28
        CPPUNIT_ASSERT_EQUAL(28.0, dEploidIO1.altCount_[20]);

        // Pf3D7_01_v3 180192 0.747258396161626
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.747258396161626,
                                     dEploidIO1.plaf_[99], 0.000000000001);
        // Pf3D7_01_v3 180192 199
        CPPUNIT_ASSERT_EQUAL(199.0, dEploidIO1.refCount_[99]);
        // Pf3D7_01_v3 180192 0
        CPPUNIT_ASSERT_EQUAL(0.0, dEploidIO1.altCount_[99]);
    }


    void testForbidMoves() {
        char *argv1[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1", "-forbidUpdateProp",
                         "-forbidUpdateSingle", "-forbidUpdatePair" };
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(14, argv1));
        DEploidIO dEploidIO1(14, argv1);
        CPPUNIT_ASSERT_EQUAL(dEploidIO1.doUpdateProp(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIO1.doUpdateSingle(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIO1.doUpdatePair(), false);

        CPPUNIT_ASSERT_NO_THROW(DEploidIO(13, argv1));
        DEploidIO dEploidIO2(13, argv1);
        CPPUNIT_ASSERT_EQUAL(dEploidIO2.doUpdateProp(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIO2.doUpdateSingle(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIO2.doUpdatePair(), true);

        CPPUNIT_ASSERT_NO_THROW(DEploidIO(12, argv1));
        DEploidIO dEploidIO3(12, argv1);
        CPPUNIT_ASSERT_EQUAL(dEploidIO3.doUpdateProp(), false);
        CPPUNIT_ASSERT_EQUAL(dEploidIO3.doUpdateSingle(), true);
        CPPUNIT_ASSERT_EQUAL(dEploidIO3.doUpdatePair(), true);

        CPPUNIT_ASSERT_NO_THROW(DEploidIO(11, argv1));
        DEploidIO dEploidIO4(11, argv1);
        CPPUNIT_ASSERT_EQUAL(dEploidIO4.doUpdateProp() , true);
        CPPUNIT_ASSERT_EQUAL(dEploidIO4.doUpdateSingle() , true);
        CPPUNIT_ASSERT_EQUAL(dEploidIO4.doUpdatePair() , true);
    }


    void testInvalidVcf() {
        char *argv1[] = { "./dEploid",
                 "-vcf", "data/testData/PG0389-C.test.vcf",
                 "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv1), InvalidInputFile);
    }


    void testVcfHeader() {
        char *argv1[] = { "./dEploid",
              "-vcf", "data/testData/crappyVcf/badHeaderFieldNames.alt.vcf",
              "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv1), VcfInvalidHeaderFieldNames);

        char *argv2[] = { "./dEploid",
              "-vcf", "data/testData/crappyVcf/badHeaderFieldNames.chrom2.vcf",
              "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv2), VcfInvalidHeaderFieldNames);

        char *argv3[] = { "./dEploid",
              "-vcf", "data/testData/crappyVcf/badHeaderFieldNames.chrom.vcf",
              "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv3), VcfInvalidHeaderFieldNames);

        char *argv4[] = { "./dEploid",
              "-vcf", "data/testData/crappyVcf/badHeaderFieldNames.filter.vcf",
              "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv4), VcfInvalidHeaderFieldNames);

        char *argv5[] = { "./dEploid",
              "-vcf", "data/testData/crappyVcf/badHeaderFieldNames.format.vcf",
              "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv5), VcfInvalidHeaderFieldNames);

        char *argv6[] = { "./dEploid",
                "-vcf", "data/testData/crappyVcf/badHeaderFieldNames.id.vcf",
                 "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv6), VcfInvalidHeaderFieldNames);

        char *argv7[] = { "./dEploid",
                 "-vcf", "data/testData/crappyVcf/badHeaderFieldNames.info.vcf",
                 "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv7), VcfInvalidHeaderFieldNames);

        char *argv8[] = { "./dEploid",
                 "-vcf", "data/testData/crappyVcf/badHeaderFieldNames.pos.vcf",
                 "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv8), VcfInvalidHeaderFieldNames);

        char *argv9[] = { "./dEploid",
                 "-vcf", "data/testData/crappyVcf/badHeaderFieldNames.qual.vcf",
                 "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv9), VcfInvalidHeaderFieldNames);

        char *argv10[] = { "./dEploid",
                 "-vcf", "data/testData/crappyVcf/badHeaderFieldNames.ref.vcf",
                 "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv10), VcfInvalidHeaderFieldNames);
    }

    void testVcfNoAD() {
        char *argv1[] = { "./dEploid",
            "-vcf", "data/testData/crappyVcf/badVariant.noAD.vcf",
            "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv1), VcfCoverageFieldNotFound);
    }

    void testInvalidVcfGz() {
        char *argv1[] = { "./dEploid",
            "-vcf", "data/testData/PG0389-C.test.vcf.gz",
            "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv1), InvalidInputFile);
    }


    void testVcfGzHeader() {
        char *argv1[] = { "./dEploid",
          "-vcf", "data/testData/crappyVcfGz/badHeaderFieldNames.alt.vcf.gz",
          "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv1), VcfInvalidHeaderFieldNames);

        char *argv2[] = { "./dEploid",
          "-vcf", "data/testData/crappyVcfGz/badHeaderFieldNames.chrom2.vcf.gz",
          "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv2), VcfInvalidHeaderFieldNames);

        char *argv3[] = { "./dEploid",
          "-vcf", "data/testData/crappyVcfGz/badHeaderFieldNames.chrom.vcf.gz",
          "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv3), VcfInvalidHeaderFieldNames);

        char *argv4[] = { "./dEploid",
          "-vcf", "data/testData/crappyVcfGz/badHeaderFieldNames.filter.vcf.gz",
          "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv4), VcfInvalidHeaderFieldNames);

        char *argv5[] = { "./dEploid",
         "-vcf", "data/testData/crappyVcfGz/badHeaderFieldNames.format.vcf.gz",
         "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv5), VcfInvalidHeaderFieldNames);

        char *argv6[] = { "./dEploid",
             "-vcf", "data/testData/crappyVcfGz/badHeaderFieldNames.id.vcf.gz",
             "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv6), VcfInvalidHeaderFieldNames);

        char *argv7[] = { "./dEploid",
            "-vcf", "data/testData/crappyVcfGz/badHeaderFieldNames.info.vcf.gz",
            "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv7), VcfInvalidHeaderFieldNames);

        char *argv8[] = { "./dEploid",
            "-vcf", "data/testData/crappyVcfGz/badHeaderFieldNames.pos.vcf.gz",
            "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv8), VcfInvalidHeaderFieldNames);

        char *argv9[] = { "./dEploid",
            "-vcf", "data/testData/crappyVcfGz/badHeaderFieldNames.qual.vcf.gz",
            "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv9), VcfInvalidHeaderFieldNames);

        char *argv10[] = { "./dEploid",
            "-vcf", "data/testData/crappyVcfGz/badHeaderFieldNames.ref.vcf.gz",
            "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv10), VcfInvalidHeaderFieldNames);
    }

    void testVcfGzNoAD() {
        char *argv1[] = { "./dEploid",
                "-vcf", "data/testData/crappyVcfGz/badVariant.noAD.vcf.gz",
                "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv1), VcfCoverageFieldNotFound);
    }

    void testVcfGzNoVQSLOD() {
        char *argv1[] = { "./dEploid",
                "-vcf", "data/testData/crappyVcf/badVariant.noVQSLOD.vcf",
                "-plaf", "data/testData/labStrains.test.PLAF.txt", "-noPanel"};
        CPPUNIT_ASSERT_THROW(DEploidIO(6, argv1), VcfVQSLODNotFound);
    }

    void testVcfOutUnSpecified() {
        char *argv1[] = { "./dEploid",
                         "-vcf", "data/testData/PG0390-C.test.vcf",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-noPanel", "-z"};
        CPPUNIT_ASSERT_THROW(DEploidIO(7, argv1), VcfOutUnSpecified);
    }

    void testChromPainting() {
        char *argv1[] = { "./dEploid",
                         "-vcf", "data/testData/PG0390-C.test.vcf",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-painting", "data/testData/PG0390-C.test.nopanel.hap",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-inbreeding"};
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(9, argv1));
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(10, argv1));
        DEploidIO dEploidIO(10, argv1);
        CPPUNIT_ASSERT_EQUAL(dEploidIO.doLsPainting(), true);
        CPPUNIT_ASSERT_NO_THROW(dEploidIO.chromPainting());

        char *argv2[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-initialP", "0.1", "0.2", "0.3", "0.4",
                         "-ibdPainting"};
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(15, argv2));
        DEploidIO dEploidIOibd(15, argv2);
        CPPUNIT_ASSERT_EQUAL(dEploidIOibd.doIbdPainting(), true);
        CPPUNIT_ASSERT_NO_THROW(dEploidIOibd.paintIBD());
    }

    void testInvalidK() {
        char *argv1[] = { "./dEploid",
                         "-vcf", "data/testData/PG0390-C.test.vcf",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-k", "1", "-ibd"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv1), InvalidK);

        char *argv2[] = { "./dEploid",
                         "-vcf", "data/testData/PG0390-C.test.vcf",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-initialP", "1", "-ibd"};
        CPPUNIT_ASSERT_THROW(DEploidIO(10, argv2), InvalidK);
    }

    void testComputeObsWsaf() {
        char *argv[] = { "./dEploid",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-lasso" };
        CPPUNIT_ASSERT_NO_THROW(DEploidIO(10, argv));
        DEploidIO dEploidIOlasso(10, argv);
        CPPUNIT_ASSERT_NO_THROW(dEploidIOlasso.dEploidLasso());
        CPPUNIT_ASSERT_EQUAL(dEploidIOlasso.lassoPanels.size(), (size_t)14);
        CPPUNIT_ASSERT_EQUAL(dEploidIOlasso.lassoPlafs.size(), (size_t)14);
        CPPUNIT_ASSERT_EQUAL(dEploidIOlasso.lassoRefCount.size(), (size_t)14);
        CPPUNIT_ASSERT_EQUAL(dEploidIOlasso.lassoAltCount.size(), (size_t)14);
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestIO);
