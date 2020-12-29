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
#include "src/utility.hpp"


class TestUtility : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE(TestUtility);
    CPPUNIT_TEST(testMinValue);
    CPPUNIT_TEST(testMaxValue);
    CPPUNIT_TEST(testVecDiff);
    CPPUNIT_TEST(testVecSum);
    CPPUNIT_TEST(testVecProd);
    CPPUNIT_TEST(testComputeCdf);
    CPPUNIT_TEST(testSampleIndexGivenProp);
    CPPUNIT_TEST(testSumOfVec);
    CPPUNIT_TEST(testSumOfMat);
    CPPUNIT_TEST(testReshapeMatToVec);
    CPPUNIT_TEST(testNormalizeBySumMat);
    CPPUNIT_TEST(testNormalPdf);
    CPPUNIT_TEST(testCalcLLK);
    CPPUNIT_TEST(testCalcLLKs);
    CPPUNIT_TEST(testBinomialPdf);
    CPPUNIT_TEST(testLogBetaPdf);
    CPPUNIT_TEST(testRBeta);

    CPPUNIT_TEST_SUITE_END();

 private:
    vector <double> vec1;
    vector <double> vec2;
    vector < vector < double> > mat1;
    MersenneTwister* rg;
    double ep1;
    double ep2;
    double ep3;
    double ep4;
    size_t nRepeat;
    double err_;
    double scalingFactor_;

    void testSampleIndexGivenPropCore(vector <double> prop) {
        vector <int> counter(prop.size(), 0);
        for (size_t i = 0; i < nRepeat; i++) {
            counter[sampleIndexGivenProp(this->rg, prop)]++;
        }
        for (size_t i = 0; i < prop.size(); i++) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(prop[i],
              (double)counter[i]/(double)nRepeat, ep1);
        }
    }

 public:
    void setUp() {
        rg = new MersenneTwister((size_t)1);
        vec1 = vector < double > ({4.0, 2.0, 3.0, 1.0});
        vec2 = vector < double > ({0.4, 0.3, 0.1, 0.2});
        mat1.push_back(vector <double> ({1.0, 2.0, 3.0}));
        mat1.push_back(vector <double> ({4.0, 5.0, 6.0}));
        mat1.push_back(vector <double> ({7.0, 8.0, 9.0}));
        ep1 = 0.001;
        ep2 = 0.0000001;
        ep4 = 0.00001;
        ep3 = 0.00000000001;
        nRepeat = 1000000;
        err_ = 0.01;
        scalingFactor_ = 100.0;
    }


    void tearDown() {
        delete rg;
    }


    void testMinValue() {
        double min1 = min_value(vec1);
        CPPUNIT_ASSERT_EQUAL(1.0, min1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1, min_value(vec2), ep3);
    }


    void testMaxValue() {
        CPPUNIT_ASSERT_EQUAL(4.0, max_value(vec1));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4, max_value(vec2), ep3);
    }


    void testNormalPdf() {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3989423,
          normal_pdf(0.0, 0.0, 1.0), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.2419707,
          normal_pdf(-1.0, 0.0, 1.0), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.2419707,
          normal_pdf(1.0, 0.0, 1.0), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0001338302,
          normal_pdf(4.0, 0.0, 1.0), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.2419707,
          normal_pdf(0.0, 1.0, 1.0), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05399097,
          normal_pdf(-1.0, 1.0, 1.0), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3989423,
          normal_pdf(1.0, 1.0, 1.0), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.004431848,
          normal_pdf(4.0, 1.0, 1.0), ep2);
    }


    void testSampleIndexGivenProp() {
        CPPUNIT_ASSERT_EQUAL((size_t)0,
          sampleIndexGivenProp(this->rg, vector <double> ({1.0, 0.0, 0.0})));
        CPPUNIT_ASSERT_EQUAL((size_t)1,
          sampleIndexGivenProp(this->rg, vector <double> ({0.0, 1.0, 0.0})));
        CPPUNIT_ASSERT_EQUAL((size_t)2,
          sampleIndexGivenProp(this->rg, vector <double> ({0.0, 0.0, 1.0})));

        this->testSampleIndexGivenPropCore(
          vector <double> ({0.9, 0.05, 0.05}));
        this->testSampleIndexGivenPropCore(
          vector <double> ({0.4, 0.6}));
        this->testSampleIndexGivenPropCore(
          vector <double> ({0.33, 0.33, 0.34}));
        this->testSampleIndexGivenPropCore(
          vector <double> ({0.1, 0.2, 0.3, 0.4}));
    }


    void testVecDiff() {
        vector <double> diff1 = vecDiff(vec1, vec2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(3.6, diff1[0], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.7, diff1[1], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2.9, diff1[2], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.8, diff1[3], ep3);

        vector <double> diff2 = vecDiff(vec2, vec1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.6, diff2[0], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.7, diff2[1], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-2.9, diff2[2], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.8, diff2[3], ep3);
    }


    void testVecSum() {
        vector <double> sum1 = vecSum(vec1, vec2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(4.4, sum1[0], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2.3, sum1[1], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(3.1, sum1[2], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.2, sum1[3], ep3);
    }


    void testVecProd() {
        vector <double> prod1 = vecProd(vec1, vec2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.6, prod1[0], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.6, prod1[1], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3, prod1[2], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.2, prod1[3], ep3);
    }


    void testComputeCdf() {
        vector <double> prop1(vec1);
        (void)normalizeBySum(prop1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4, prop1[0], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.2, prop1[1], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3, prop1[2], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1, prop1[3], ep3);

        vector <double> cdf1 = computeCdf(prop1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4, cdf1[0], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.6, cdf1[1], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.9, cdf1[2], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, cdf1[3], ep3);

        vector <double> prop2(vec2);
        (void)normalizeBySum(prop2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4, prop2[0], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3, prop2[1], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1, prop2[2], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.2, prop2[3], ep3);

        vector <double> cdf2 = computeCdf(prop2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4, cdf2[0], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7, cdf2[1], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.8, cdf2[2], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, cdf2[3], ep3);
    }


    void testSumOfMat() {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(45.0, sumOfMat(mat1), ep3);
    }


    void testSumOfVec() {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, sumOfVec(vec1), ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, sumOfVec(vec2), ep3);
    }


    void testNormalizeBySumMat() {
        vector < vector <double>> tmpMat1 = mat1;
        (void)normalizeBySumMat(tmpMat1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02222222, tmpMat1[0][0], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.04444444, tmpMat1[0][1], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.06666667, tmpMat1[0][2], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08888889, tmpMat1[1][0], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.11111111, tmpMat1[1][1], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.13333333, tmpMat1[1][2], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.15555556, tmpMat1[2][0], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.17777778, tmpMat1[2][1], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.20000000, tmpMat1[2][2], ep2);
    }


    void testReshapeMatToVec() {
        vector <double> tmpvec1 = reshapeMatToVec(mat1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tmpvec1[0], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, tmpvec1[1], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, tmpvec1[2], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, tmpvec1[3], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, tmpvec1[4], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, tmpvec1[5], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(7.0, tmpvec1[6], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(8.0, tmpvec1[7], ep3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0, tmpvec1[8], ep3);
    }


    void testCalcLLK() {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0,
          calcLLK(0, 0, 0.0, this->err_, this->scalingFactor_), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0,
          calcLLK(0, 0, 1.0, this->err_, this->scalingFactor_), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0,
          calcLLK(0, 0, 0.5, this->err_, this->scalingFactor_), ep2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.09622803,
          calcLLK(10, 0, 0.0, this->err_, this->scalingFactor_), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-31.38367782,
          calcLLK(10, 0, 1.0, this->err_, this->scalingFactor_), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-6.520005876,
          calcLLK(10, 0, 0.5, this->err_, this->scalingFactor_), ep2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(-92.39749823,
          calcLLK(0, 50, 0.0, this->err_, this->scalingFactor_), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.4088264558,
          calcLLK(0, 50, 1.0, this->err_, this->scalingFactor_), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-26.30680376,
          calcLLK(0, 50, 0.5, this->err_, this->scalingFactor_), ep2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(-39.454802987,
          calcLLK(11, 13, 0.0, this->err_, this->scalingFactor_), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-82.7561739077,
          calcLLK(41, 2, 1.0, this->err_, this->scalingFactor_), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-71.6257227,
          calcLLK(23, 99, 0.5, this->err_, this->scalingFactor_), ep2);
    }


    void testCalcLLKs() {
        vector <double> ref({ 0, 0, 0, 10, 10, 10, 0, 0, 0, 11, 41, 23});
        vector <double> alt({ 0, 0, 0, 0, 0, 0, 50, 50, 50, 13, 2, 99});
        vector <double> wsaf({
          .0, 1.0, .5, .0, 1.0, .5, .0, 1.0, .5, .0, 1.0, .5});

        vector <double> llk1 =
          calcLLKs(ref, alt, wsaf, 0, 3, this->scalingFactor_);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, llk1[0], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, llk1[0], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, llk1[0], ep2);

        vector <double> llk2 =
          calcLLKs(ref, alt, wsaf, 6, 3, this->scalingFactor_);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-92.39749823, llk2[0], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.4088264558, llk2[1], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-26.30680376, llk2[2], ep2);

        vector <double> llk3 =
          calcLLKs(ref, alt, wsaf, 3, 7, this->scalingFactor_);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.09622803, llk3[0], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-31.38367782, llk3[1], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-6.520005876, llk3[2], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-92.39749823, llk3[3], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.4088264558, llk3[4], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-26.30680376, llk3[5], ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-39.454802987, llk3[6], ep2);
    }


    void testLogBetaPdf() {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.442989, betaPdf(.1, .1, .1), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.8142104, logBetaPdf(.1, .1, .1), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.789602, betaPdf(.1, .1, .9), ep2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.2362263, logBetaPdf(.1, .1, .9), ep2);
    }


    void testBinomialPdf() {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4444444,
          binomialPdf(0, 2, (1.0/3.0)), ep4);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4444444,
          binomialPdf(1, 2, (1.0/3.0)), ep4);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.111111,
          binomialPdf(static_cast<int>(2),
            static_cast<int>(2), (1.0/3.0)), ep4);
    }


    void testRBetaCore(double a, double b) {
        double Ex = 0;
        double Ex2 = 0;
        double tmp;
        for (size_t i = 0; i < this->nRepeat; i++) {
            tmp = rBeta(a, b, rg);
            Ex += tmp;
            Ex2 += pow(tmp, 2.0);
        }
        Ex /= nRepeat;
        Ex2 /= nRepeat;
        double mean = a/(a+b);
        double var = (a*b) / (pow((a+b), 2.0) * (a+b+1));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(mean, Ex, ep1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(var, Ex2-pow(Ex, 2.0), ep1);
    }

    void testRBeta() {
        this->testRBetaCore(4.0, 3.0);
        this->testRBetaCore(14.0, 31.0);
        this->testRBetaCore(53.0, 63.0);
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestUtility);
