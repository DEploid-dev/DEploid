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

class TestWorkflow : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE(TestWorkflow);
    CPPUNIT_TEST(testWorkflowBest);
    CPPUNIT_TEST_SUITE_END();

 public:
    void setUp() {
    }

    void tearDown() {
    }

    void testWorkflowBest(){
        char *argv1[] = { "./dEploid",
                         "-vcf", "data/testData/PG0390-C.test.vcf.gz",
                         "-sample", "PG0390-C", "-plafFromVcf",
                         "-exclude", "data/testData/labStrains.test.exclude.txt.gz"
                         "-panel", "data/testData/labStrains.test.panel.txt.gz",
                         "-best"};
        DEploidIO tmp(11, argv1);
        CPPUNIT_ASSERT_NO_THROW(tmp.workflow_best());

    }
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestWorkflow);

