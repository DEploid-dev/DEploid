#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "updateHap.hpp"


class TestUpdateHap : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestUpdateHap );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testVirtualFunctions );
    CPPUNIT_TEST_SUITE_END();

  private:
    UpdateHap * updateHapPlaf_;
    UpdateHap * updateHapPanel_;

    vector <double> refCount_;
    vector <double> altCount_;
    vector <double> plaf_;
    vector <double> expectedWsaf_;
    vector <double> proportion_;
    vector < vector <double> > haplotypes_;
    MersenneTwister* rg_;
    size_t segmentStartIndex_;
    size_t nLoci_;
    Panel* panel_;

  public:
    void setUp(){
        this->rg_ = new MersenneTwister((size_t)1);
        this->panel_ = new Panel();

        this->updateHapPlaf_ = new UpdateHap( refCount_,
                                          altCount_,
                                          plaf_,
                                          expectedWsaf_,
                                          proportion_,
                                          haplotypes_,
                                          rg_,
                                          segmentStartIndex_,
                                          nLoci_,
                                          NULL );

        this->updateHapPanel_ = new UpdateHap( refCount_,
                                          altCount_,
                                          plaf_,
                                          expectedWsaf_,
                                          proportion_,
                                          haplotypes_,
                                          rg_,
                                          segmentStartIndex_,
                                          nLoci_,
                                          panel_ );
    }


    void tearDown(){
        delete updateHapPlaf_;
        delete updateHapPanel_;
        delete panel_;
        this->rg_->clearFastFunc();
        delete rg_;
    }


    void testMainConstructor(){ }


    void testVirtualFunctions(){
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->calcExpectedWsaf(expectedWsaf_, proportion_, haplotypes_ ) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->calcHapLLKs(refCount_, altCount_) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->buildEmission(0.1) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->samplePaths() );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->addMissCopying(0.1) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->updateLLK() );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->sampleHapIndependently(plaf_) );

        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPanel_->calcExpectedWsaf(expectedWsaf_, proportion_, haplotypes_ ) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPanel_->calcHapLLKs(refCount_, altCount_) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPanel_->buildEmission(0.1) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPanel_->samplePaths() );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPanel_->addMissCopying(0.1) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPanel_->updateLLK() );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPanel_->sampleHapIndependently(plaf_) );
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestUpdateHap );


class TestUpdateSingleHap : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestUpdateSingleHap );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST ( testEmissionProb0 );
    CPPUNIT_TEST ( testEmissionBasicVersion0 );
    CPPUNIT_TEST_SUITE_END();

  private:
    UpdateSingleHap * updateSingleHapPanel_;
    UpdateSingleHap * updateSingleHapPlaf_;
    vector <double> refCount_;
    vector <double> altCount_;
    vector <double> plaf_;
    vector <double> expectedWsaf_;
    vector <double> proportion_;
    vector < vector <double> > haplotypes_;
    MersenneTwister* rg_;
    size_t segmentStartIndex_;
    size_t nLoci_;
    Panel* panel_;
    double missCopyProb_;
    size_t strainIndex_;

  public:

    void setUp(){
        this->expectedWsaf_ = vector <double> ();
        this->rg_ = new MersenneTwister((size_t)1);
        this->panel_ = new Panel();
        this->nLoci_ = 0;

        this->updateSingleHapPlaf_ = new UpdateSingleHap( refCount_,
                                          altCount_,
                                          plaf_,
                                          expectedWsaf_,
                                          proportion_,
                                          haplotypes_,
                                          rg_,
                                          segmentStartIndex_,
                                          nLoci_,
                                          NULL,
                                          missCopyProb_,
                                          strainIndex_ );

        this->updateSingleHapPanel_ = new UpdateSingleHap( refCount_,
                                          altCount_,
                                          plaf_,
                                          expectedWsaf_,
                                          proportion_,
                                          haplotypes_,
                                          rg_,
                                          segmentStartIndex_,
                                          nLoci_,
                                          panel_,
                                          missCopyProb_,
                                          strainIndex_ );
}


    void tearDown(){
        delete updateSingleHapPlaf_;
        delete updateSingleHapPanel_;
        delete panel_;
        this->rg_->clearFastFunc();
        delete rg_;
    }

    void testMainConstructor(){ }

    void testEmissionProb0 (){
        double llk0_1 = log(0.05), llk0_2 = log(0.03), llk0_3 = log(0);
        double llk1_1 = log(0.5), llk1_2 = log(0.0003), llk1_3 = log(1);
        this->updateSingleHapPanel_->llk0_ = vector <double> ({llk0_1, llk0_2, llk0_3});
        this->updateSingleHapPanel_->llk1_ = vector <double> ({llk1_1, llk1_2, llk1_3});
        this->updateSingleHapPanel_->nLoci_ = 3;
        this->updateSingleHapPanel_->buildEmission ( 0.0 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05/0.5, updateSingleHapPanel_->emission_[0][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel_->emission_[0][1], 0.000000000001);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel_->emission_[1][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0003/0.03, updateSingleHapPanel_->emission_[1][1], 0.000000000001);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, updateSingleHapPanel_->emission_[2][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel_->emission_[2][1], 0.000000000001);
    }


    void testEmissionBasicVersion0 (){
        double llk0_1 = log(0.05), llk0_2 = log(0.03), llk0_3 = log(0);
        double llk1_1 = log(0.5), llk1_2 = log(0.0003), llk1_3 = log(1);
        this->updateSingleHapPanel_->llk0_ = vector <double> ({llk0_1, llk0_2, llk0_3});
        this->updateSingleHapPanel_->llk1_ = vector <double> ({llk1_1, llk1_2, llk1_3});
        this->updateSingleHapPanel_->nLoci_ = 3;
        this->updateSingleHapPanel_->buildEmissionBasicVersion ( 0.0 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_1), updateSingleHapPanel_->emission_[0][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_1), updateSingleHapPanel_->emission_[0][1], 0.000000000001);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_2), updateSingleHapPanel_->emission_[1][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_2), updateSingleHapPanel_->emission_[1][1], 0.000000000001);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_3), updateSingleHapPanel_->emission_[2][0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_3), updateSingleHapPanel_->emission_[2][1], 0.000000000001);
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestUpdateSingleHap );
