#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "ibd.hpp"

class TestIBDUtility: public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestIBDUtility );
    CPPUNIT_TEST( testConvertIntToBinary );
    CPPUNIT_TEST( testNchoose2 );
    CPPUNIT_TEST_SUITE_END();

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


class TestIBDconfig : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestIBDconfig );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testCountEffectiveK );
    CPPUNIT_TEST_SUITE_END();

  private:
    IBDconfiguration* ibd3_;
    IBDconfiguration* ibd4_;
    IBDconfiguration* ibd5_;

  public:
    void setUp() {
        ibd3_ = new IBDconfiguration();
        CPPUNIT_ASSERT_NO_THROW( ibd3_->buildIBDconfiguration(3));
        ibd4_ = new IBDconfiguration();
        CPPUNIT_ASSERT_NO_THROW( ibd4_->buildIBDconfiguration(4));
        ibd5_ = new IBDconfiguration();
        CPPUNIT_ASSERT_NO_THROW( ibd5_->buildIBDconfiguration(5));
    }


    void tearDown() {
        delete ibd3_;
        delete ibd4_;
        delete ibd5_;
    }


    vector <int> computeEffectiveKCount(int k, IBDconfiguration* ibd){
        vector <int> ret(k);
        for ( int index = 0; index < k; index++){
            for (int effk : ibd->effectiveK ){
                if (effk == (index+1)){
                    ret[index]++;
                }
            }
        }
        return ret;
    }



    void testMainConstructor(){
        CPPUNIT_ASSERT_NO_THROW( IBDconfiguration());
        //CPPUNIT_ASSERT_NO_THROW( IBDconfiguration(5));
        CPPUNIT_ASSERT_EQUAL( ibd3_->states.size(), (size_t)5);
        CPPUNIT_ASSERT_EQUAL( ibd4_->states.size(), (size_t)15);
        CPPUNIT_ASSERT_EQUAL( ibd5_->states.size(), (size_t)52);
    }


    void testCountEffectiveK(){
        // k = 3
        vector <int> count3 = computeEffectiveKCount(3, ibd3_);
        //table(make.ibd.mat.joe(3)$k.eff)
        //1 2 3
        //1 3 1
        CPPUNIT_ASSERT_EQUAL( count3[0], (int)1); // 1
        CPPUNIT_ASSERT_EQUAL( count3[1], (int)3); // 2
        CPPUNIT_ASSERT_EQUAL( count3[2], (int)1); // 3

        // k = 4
        vector <int> count4 = computeEffectiveKCount(4, ibd4_);
        //table(make.ibd.mat.joe(4)$k.eff)
        //1 2 3 4
        //1 7 6 1
        CPPUNIT_ASSERT_EQUAL( count4[0], (int)1); // 1
        CPPUNIT_ASSERT_EQUAL( count4[1], (int)7); // 2
        CPPUNIT_ASSERT_EQUAL( count4[2], (int)6); // 3
        CPPUNIT_ASSERT_EQUAL( count4[3], (int)1); // 4

        // k = 5
        vector <int> count5 = computeEffectiveKCount(5, ibd5_);
        //table(make.ibd.mat.joe(5)$k.eff)
        //1  2  3  4  5
        //1 15 25 10  1
        CPPUNIT_ASSERT_EQUAL( count5[0], (int)1);  // 1
        CPPUNIT_ASSERT_EQUAL( count5[1], (int)15); // 2
        CPPUNIT_ASSERT_EQUAL( count5[2], (int)25); // 3
        CPPUNIT_ASSERT_EQUAL( count5[3], (int)10); // 4
        CPPUNIT_ASSERT_EQUAL( count5[4], (int)1);  // 5
    }


};

class TestHprior : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestHprior );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST_SUITE_END();

  private:
    Hprior* hprior3_;
    Hprior* hprior4_;
    Hprior* hprior5_;
    vector <double> plaf;
  public:
    void setUp() {
        plaf = vector <double> ({0.0, 0.5, 1.0});
        hprior3_ = new Hprior;
        CPPUNIT_ASSERT_NO_THROW(hprior3_->buildHprior(3, plaf));
        hprior4_ = new Hprior;(4, plaf);
        CPPUNIT_ASSERT_NO_THROW(hprior4_->buildHprior(4, plaf));
        hprior5_ = new Hprior;(5, plaf);
        CPPUNIT_ASSERT_NO_THROW(hprior5_->buildHprior(5, plaf));
    }


    void tearDown() {
        delete hprior3_;
        delete hprior4_;
        delete hprior5_;
    }


    void testMainConstructor(){
        CPPUNIT_ASSERT_EQUAL( hprior3_->nState(), (size_t)22); // This is equal to sum(2^IBDconfiguration(3).effectiveK)
        CPPUNIT_ASSERT_EQUAL( hprior4_->nState(), (size_t)94);
        CPPUNIT_ASSERT_EQUAL( hprior5_->nState(), (size_t)454);
    }


};

CPPUNIT_TEST_SUITE_REGISTRATION(TestIBDUtility);
CPPUNIT_TEST_SUITE_REGISTRATION(TestIBDconfig);
CPPUNIT_TEST_SUITE_REGISTRATION(TestHprior);
