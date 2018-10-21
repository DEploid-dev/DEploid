/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
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
#include "src/vcfReader.hpp"

class TestVCF : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE(TestVCF);
    CPPUNIT_TEST(testMainConstructor);
    CPPUNIT_TEST_SUITE_END();

 private:
    VcfReader* vcf_;
    VcfReader* vcfGz_;
    double eps;

 public:
    void setUp() {
        this->vcf_ = new VcfReader("data/testData/PG0390-C.test.vcf");
        this->vcfGz_ = new VcfReader("data/testData/PG0390-C.test.vcf.gz");
        this->eps = 0.00000000001;
    }

    void tearDown() {
        delete this->vcf_;
        delete this->vcfGz_;
    }

    void testMainConstructor() {
        CPPUNIT_ASSERT_NO_THROW(this->vcf_->finalize());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(8.08, this->vcf_->vqslod[0], this->eps);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.617, this->vcf_->vqslod[1], this->eps);
        CPPUNIT_ASSERT_EQUAL(this->vcf_->vqslod.size(),
                             this->vcf_->refCount.size());
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestVCF);
