#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "utility.hpp"


class TestUtility : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestUtility );
    CPPUNIT_TEST( testMinValue );
    CPPUNIT_TEST( testMaxValue );
    CPPUNIT_TEST( testVecDiff );
    CPPUNIT_TEST( testVecSum );
    CPPUNIT_TEST( testVecProd );
    CPPUNIT_TEST( testComputeCdf );
    CPPUNIT_TEST( testSampleIndexGivenProp );
    CPPUNIT_TEST( testSumOfVec );
    CPPUNIT_TEST( testSumOfMat );
    CPPUNIT_TEST( testReshapeMatToVec );
    CPPUNIT_TEST( testNormalizeBySumMat );
    CPPUNIT_TEST( testNormalPdf );
    CPPUNIT_TEST( testCalcLLK );
    CPPUNIT_TEST( testCalcLLKs );
    CPPUNIT_TEST( testConvertIntToBinary );
    CPPUNIT_TEST( testNchoose2 );
    CPPUNIT_TEST_SUITE_END();

  private:
    vector <double> vec1;
    vector <double> vec2;
    vector < vector < double> > mat1;
    MersenneTwister* rg;
    double epsilon1;
    double epsilon2;
    double epsilon3;
    size_t nRepeat;

    void testSampleIndexGivenPropCore( vector <double> prop ){
        vector <int> counter ( prop.size(), 0 );
        for ( size_t i = 0; i < nRepeat; i++ ){
            counter[sampleIndexGivenProp( this->rg, prop)]++;
        }
        for ( size_t i = 0; i < prop.size(); i++ ){
            CPPUNIT_ASSERT_DOUBLES_EQUAL ( prop[i], (double)counter[i]/(double)nRepeat, epsilon1 );
        }
    }

  public:
    void setUp() {
        rg = new MersenneTwister((size_t)1);
        vec1 = vector < double > ({4.0, 2.0, 3.0, 1.0});
        vec2 = vector < double > ({0.4, 0.3, 0.1, 0.2});
        mat1.push_back( vector <double> ({1.0, 2.0, 3.0}) );
        mat1.push_back( vector <double> ({4.0, 5.0, 6.0}) );
        mat1.push_back( vector <double> ({7.0, 8.0, 9.0}) );
        epsilon1 = 0.001;
        epsilon2 = 0.0000001;
        epsilon3 = 0.00000000001;
        nRepeat = 1000000;
    }


    void tearDown() {
        delete rg;
    }


    void testMinValue(){
        double min1 =  min_value(vec1);
        CPPUNIT_ASSERT_EQUAL ( 1.0, min1 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.1, min_value(vec2), epsilon3 );
    }


    void testMaxValue(){
        CPPUNIT_ASSERT_EQUAL ( 4.0, max_value(vec1) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.4, max_value(vec2), epsilon3 );
    }


    void testNormalPdf(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3989423, normal_pdf(0.0, 0.0, 1.0), epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.2419707, normal_pdf(-1.0, 0.0, 1.0), epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.2419707, normal_pdf(1.0, 0.0, 1.0), epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.0001338302, normal_pdf(4.0, 0.0, 1.0), epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.2419707, normal_pdf(0.0, 1.0, 1.0), epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.05399097, normal_pdf(-1.0, 1.0, 1.0), epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3989423, normal_pdf(1.0, 1.0, 1.0), epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.004431848, normal_pdf(4.0, 1.0, 1.0), epsilon2 );
    }


    void testSampleIndexGivenProp(){
        CPPUNIT_ASSERT_EQUAL ( (size_t)0, sampleIndexGivenProp( this->rg, vector <double> ({1.0, 0.0, 0.0})) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)1, sampleIndexGivenProp( this->rg, vector <double> ({0.0, 1.0, 0.0})) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)2, sampleIndexGivenProp( this->rg, vector <double> ({0.0, 0.0, 1.0})) );

        this->testSampleIndexGivenPropCore( vector <double> ({0.9, 0.05, 0.05}) );
        this->testSampleIndexGivenPropCore( vector <double> ({0.4, 0.6}) );
        this->testSampleIndexGivenPropCore( vector <double> ({0.33, 0.33, 0.34}) );
        this->testSampleIndexGivenPropCore( vector <double> ({0.1, 0.2, 0.3, 0.4}) );
    }


    void testVecDiff(){
        vector <double> diff1 = vecDiff(vec1, vec2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 3.6, diff1[0], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1.7, diff1[1], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.9, diff1[2], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.8, diff1[3], epsilon3 );

        vector <double> diff2 = vecDiff(vec2, vec1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( -3.6, diff2[0], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( -1.7, diff2[1], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( -2.9, diff2[2], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( -0.8, diff2[3], epsilon3 );
    }


    void testVecSum(){
        vector <double> sum1 = vecSum(vec1, vec2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 4.4, sum1[0], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.3, sum1[1], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 3.1, sum1[2], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1.2, sum1[3], epsilon3 );
    }


    void testVecProd(){
        vector <double> prod1 = vecProd(vec1, vec2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1.6, prod1[0], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.6, prod1[1], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, prod1[2], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.2, prod1[3], epsilon3 );
    }


    void testComputeCdf(){
        vector <double> prop1 (vec1);
        (void)normalizeBySum(prop1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.4, prop1[0], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.2, prop1[1], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, prop1[2], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.1, prop1[3], epsilon3 );

        vector <double> cdf1 = computeCdf(prop1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.4, cdf1[0], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.6, cdf1[1], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.9, cdf1[2], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1.0, cdf1[3], epsilon3 );

        vector <double> prop2 (vec2);
        (void)normalizeBySum(prop2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.4, prop2[0], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, prop2[1], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.1, prop2[2], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.2, prop2[3], epsilon3 );

        vector <double> cdf2 = computeCdf(prop2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.4, cdf2[0], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.7, cdf2[1], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.8, cdf2[2], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1.0, cdf2[3], epsilon3 );
    }


    void testSumOfMat(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(45.0, sumOfMat(mat1), epsilon3 );
    }


    void testSumOfVec(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, sumOfVec(vec1), epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, sumOfVec(vec2), epsilon3 );
    }


    void testNormalizeBySumMat(){
        vector < vector <double>> tmpMat1 = mat1;
        (void)normalizeBySumMat (tmpMat1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.02222222, tmpMat1[0][0], epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.04444444, tmpMat1[0][1], epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.06666667, tmpMat1[0][2], epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.08888889, tmpMat1[1][0], epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.11111111, tmpMat1[1][1], epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.13333333, tmpMat1[1][2], epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.15555556, tmpMat1[2][0], epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.17777778, tmpMat1[2][1], epsilon2 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.20000000, tmpMat1[2][2], epsilon2 );
    }


    void testReshapeMatToVec(){
        vector <double> tmpvec1 = reshapeMatToVec(mat1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1.0, tmpvec1[0], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.0, tmpvec1[1], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 3.0, tmpvec1[2], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 4.0, tmpvec1[3], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 5.0, tmpvec1[4], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 6.0, tmpvec1[5], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 7.0, tmpvec1[6], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 8.0, tmpvec1[7], epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 9.0, tmpvec1[8], epsilon3 );
    }


    void testCalcLLK(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, calcLLK(0, 0, 0.0), epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, calcLLK(0, 0, 1.0), epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, calcLLK(0, 0, 0.5), epsilon2) ;

        CPPUNIT_ASSERT_DOUBLES_EQUAL (-0.09622803, calcLLK(10, 0, 0.0), epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-31.38367782, calcLLK(10, 0, 1.0), epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-6.520005876, calcLLK(10, 0, 0.5), epsilon2) ;

        CPPUNIT_ASSERT_DOUBLES_EQUAL (-92.39749823, calcLLK(0, 50, 0.0), epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-0.4088264558, calcLLK(0, 50, 1.0), epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-26.30680376, calcLLK(0, 50, 0.5), epsilon2) ;

        CPPUNIT_ASSERT_DOUBLES_EQUAL (-39.454802987, calcLLK(11, 13, 0.0), epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-82.7561739077, calcLLK(41, 2, 1.0), epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-71.6257227, calcLLK(23, 99, 0.5), epsilon2) ;
    }


    void testCalcLLKs(){
        vector <double> ref ({ 0, 0, 0, 10, 10, 10, 0, 0, 0, 11, 41, 23});
        vector <double> alt ({ 0, 0, 0, 0, 0, 0, 50, 50, 50, 13, 2, 99});
        vector <double> wsaf ({.0, 1.0, .5, .0, 1.0, .5, .0, 1.0, .5, .0, 1.0, .5 });

        vector <double> llk1 = calcLLKs( ref, alt, wsaf, 0, 3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, llk1[0], epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, llk1[0], epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, llk1[0], epsilon2) ;

        vector <double> llk2 = calcLLKs( ref, alt, wsaf, 6, 3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-92.39749823, llk2[0], epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-0.4088264558, llk2[1], epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-26.30680376, llk2[2], epsilon2) ;

        vector <double> llk3 = calcLLKs( ref, alt, wsaf, 3, 7 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-0.09622803, llk3[0], epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-31.38367782, llk3[1], epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-6.520005876, llk3[2], epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-92.39749823, llk3[3], epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-0.4088264558, llk3[4], epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-26.30680376, llk3[5], epsilon2) ;
        CPPUNIT_ASSERT_DOUBLES_EQUAL (-39.454802987, llk3[6], epsilon2) ;
    }

    void testConvertIntToBinary(){
        CPPUNIT_ASSERT_THROW (convertIntToBinary(20, 3), OutOfVectorSize);
        CPPUNIT_ASSERT_THROW (convertIntToBinary(15, 3), OutOfVectorSize);
        CPPUNIT_ASSERT_THROW (convertIntToBinary(8, 3), OutOfVectorSize);
        CPPUNIT_ASSERT_THROW_MESSAGE( "Out of vector size!", convertIntToBinary(8, 3), OutOfVectorSize);

        CPPUNIT_ASSERT_NO_THROW (convertIntToBinary(7, 3));

        vector <int> test1 = convertIntToBinary(3, 4);
        CPPUNIT_ASSERT_EQUAL (test1.size(), (size_t)4);
        CPPUNIT_ASSERT_EQUAL (test1[0], (int)0);
        CPPUNIT_ASSERT_EQUAL (test1[1], (int)0);
        CPPUNIT_ASSERT_EQUAL (test1[2], (int)1);
        CPPUNIT_ASSERT_EQUAL (test1[3], (int)1);

        CPPUNIT_ASSERT_NO_THROW(enumerateBinaryMatrixOfK(10));
        vector < vector <int> > test2 = enumerateBinaryMatrixOfK(4);
        CPPUNIT_ASSERT_EQUAL (test2.size(), (size_t)16);
        for (size_t i = 0; i < 16; i++){
            CPPUNIT_ASSERT_EQUAL (test2[i].size(), (size_t)4);
            CPPUNIT_ASSERT_EQUAL ((int)(8*test2[i][0] +
                                        4*test2[i][1] +
                                        2*test2[i][2] +
                                        1*test2[i][3]), (int)i);
        }
    }

    void testNchoose2(){
        CPPUNIT_ASSERT_THROW (nchoose2(-1), InvalidInput);
        CPPUNIT_ASSERT_THROW (nchoose2(0), InvalidInput);
        CPPUNIT_ASSERT_THROW (nchoose2(1), InvalidInput);
        CPPUNIT_ASSERT_NO_THROW (nchoose2(2));
        CPPUNIT_ASSERT_EQUAL (nchoose2(2), (int)1);
        CPPUNIT_ASSERT_EQUAL (nchoose2(3), (int)3);
        CPPUNIT_ASSERT_EQUAL (nchoose2(4), (int)6);
        CPPUNIT_ASSERT_EQUAL (nchoose2(5), (int)10);
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestUtility );
