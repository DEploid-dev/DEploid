#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "updateHap.hpp"


class TestUpdateHap : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestUpdateHap );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testVirtualFunctions );
    CPPUNIT_TEST_SUITE_END();

  private:
    UpdateHap * updateHap_;
    Panel* panel_;
    MersenneTwister* rg_;
    vector <double> refCount;
    vector <double> altCount;
    vector <double> plaf;
    vector <double> expectedWsaf;
    vector <double> proportion;
    vector < vector <double> > haplotypes;
    size_t segmentStartIndex;
    size_t nLoci;

  public:
    void setUp(){
        this->rg_ = new MersenneTwister((size_t)1);
        this->panel_ = new Panel();

        this->updateHap_ = new UpdateHap( refCount,
                                          altCount,
                                          plaf,
                                          expectedWsaf,
                                          proportion,
                                          haplotypes,
                                          rg_,
                                          segmentStartIndex,
                                          nLoci,
                                          panel_ );

    }


    void tearDown(){
        delete updateHap_;
        delete panel_;
        this->rg_->clearFastFunc();
        delete rg_;
    }


    void testMainConstructor(){
        UpdateHap tmp();
    }


    void testVirtualFunctions(){
        CPPUNIT_ASSERT_NO_THROW ( this->updateHap_->calcExpectedWsaf(expectedWsaf, proportion, haplotypes ) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHap_->calcHapLLKs(refCount, altCount) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHap_->buildEmission(0.1) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHap_->samplePaths() );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHap_->addMissCopying(0.1) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHap_->updateLLK() );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHap_->sampleHapIndependently(plaf) );
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestUpdateHap );


class TestUpdateSingleHap : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestUpdateSingleHap );
    CPPUNIT_TEST ( testEmissionProb0 );
    CPPUNIT_TEST ( testEmissionBasicVersion0 );
    CPPUNIT_TEST_SUITE_END();

  private:
    UpdateSingleHap * updateSingleHap_;

  public:
    void setUp() {
        this->updateSingleHap_ = new UpdateSingleHap();
    }

    void tearDown() {
        delete updateSingleHap_;
    }

    void testEmissionProb0 (){
        this->updateSingleHap_->emission_.clear();

        double llk0_1 = log(0.05), llk0_2 = log(0.03), llk0_3 = log(0);
        double llk1_1 = log(0.5), llk1_2 = log(0.0003), llk1_3 = log(1);
        this->updateSingleHap_->llk0_ = vector <double> ({llk0_1, llk0_2, llk0_3});
        this->updateSingleHap_->llk1_ = vector <double> ({llk1_1, llk1_2, llk1_3});
        this->updateSingleHap_->nLoci_ = 3;
        this->updateSingleHap_->buildEmission ( 0.0 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05/0.5, updateSingleHap_->emission_[0][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHap_->emission_[0][1], 0.000000000001);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHap_->emission_[1][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0003/0.03, updateSingleHap_->emission_[1][1], 0.000000000001);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, updateSingleHap_->emission_[2][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHap_->emission_[2][1], 0.000000000001);
    }


    void testEmissionBasicVersion0 (){
        this->updateSingleHap_->emission_.clear();

        double llk0_1 = log(0.05), llk0_2 = log(0.03), llk0_3 = log(0);
        double llk1_1 = log(0.5), llk1_2 = log(0.0003), llk1_3 = log(1);
        this->updateSingleHap_->llk0_ = vector <double> ({llk0_1, llk0_2, llk0_3});
        this->updateSingleHap_->llk1_ = vector <double> ({llk1_1, llk1_2, llk1_3});
        this->updateSingleHap_->nLoci_ = 3;
        this->updateSingleHap_->buildEmissionBasicVersion ( 0.0 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_1), updateSingleHap_->emission_[0][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_1), updateSingleHap_->emission_[0][1], 0.000000000001);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_2), updateSingleHap_->emission_[1][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_2), updateSingleHap_->emission_[1][1], 0.000000000001);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_3), updateSingleHap_->emission_[2][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_3), updateSingleHap_->emission_[2][1], 0.000000000001);
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestUpdateSingleHap );
