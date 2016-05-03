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
    double missCopyProb_;

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
                                          NULL, missCopyProb_ );

        this->updateHapPanel_ = new UpdateHap( refCount_,
                                          altCount_,
                                          plaf_,
                                          expectedWsaf_,
                                          proportion_,
                                          haplotypes_,
                                          rg_,
                                          segmentStartIndex_,
                                          nLoci_,
                                          panel_, missCopyProb_ );
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
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->core ( refCount_, altCount_, plaf_, expectedWsaf_, proportion_, haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->calcExpectedWsaf(expectedWsaf_, proportion_, haplotypes_ ) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->calcHapLLKs(refCount_, altCount_) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->buildEmission(0.1) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->samplePaths() );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->addMissCopying(0.1) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->updateLLK() );
        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPlaf_->sampleHapIndependently(plaf_) );

        CPPUNIT_ASSERT_NO_THROW ( this->updateHapPanel_->core ( refCount_, altCount_, plaf_, expectedWsaf_, proportion_, haplotypes_ ) );
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
    CPPUNIT_TEST ( testUpdateLLK );
    CPPUNIT_TEST ( testExpectedWsaf );
    CPPUNIT_TEST ( testCalcHapLLKs );
    CPPUNIT_TEST ( testSampleHapIndependently );
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
    double epsilon3;

  public:

    void setUp(){
        epsilon3 = 0.00000000001;

        this->refCount_ = vector <double> (  { 100, 10 , 50 , 30 , 100, 7, 50 } );
        this->altCount_ = vector <double> (  { 2,  100 , 50 , 30 , 5,  70, 30 } );

        this->haplotypes_.push_back( vector <double> ({1.0,    1.0,    1.0,    1.0,    0.0}) );
        this->haplotypes_.push_back( vector <double> ({0.0,    1.0,    0.0,    1.0,    0.0}) );
        this->haplotypes_.push_back( vector <double> ({1.0,    0.0,    0.0,    1.0,    1.0}) );
        this->haplotypes_.push_back( vector <double> ({1.0,    0.0,    0.0,    0.0,    1.0}) );
        this->haplotypes_.push_back( vector <double> ({1.0,    0.0,    1.0,    0.0,    0.0}) );
        this->haplotypes_.push_back( vector <double> ({0.0,    1.0,    1.0,    0.0,    1.0}) );
        this->haplotypes_.push_back( vector <double> ({1.0,    0.0,    1.0,    0.0,    0.0}) );

        this->proportion_ = vector <double> ({0.1, .2, .2, .15, .35} );
        this->expectedWsaf_ = vector <double> ({.65, .35, .6, .45, .3, .75, .3});

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
        CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPanel_->buildEmission ( 0.0 ) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_1-llk1_1), updateSingleHapPanel_->emission_[0][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05/0.5, updateSingleHapPanel_->emission_[0][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel_->emission_[0][1], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel_->emission_[1][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0003/0.03, updateSingleHapPanel_->emission_[1][1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_2-llk0_2), updateSingleHapPanel_->emission_[1][1], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_3-llk1_3), updateSingleHapPanel_->emission_[2][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, updateSingleHapPanel_->emission_[2][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel_->emission_[2][1], epsilon3);
    }


    void testEmissionBasicVersion0 (){
        double llk0_1 = log(0.05), llk0_2 = log(0.03), llk0_3 = log(0);
        double llk1_1 = log(0.5), llk1_2 = log(0.0003), llk1_3 = log(1);
        this->updateSingleHapPanel_->llk0_ = vector <double> ({llk0_1, llk0_2, llk0_3});
        this->updateSingleHapPanel_->llk1_ = vector <double> ({llk1_1, llk1_2, llk1_3});
        this->updateSingleHapPanel_->nLoci_ = 3;
        CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPanel_->buildEmissionBasicVersion ( 0.0 ) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_1), updateSingleHapPanel_->emission_[0][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_1), updateSingleHapPanel_->emission_[0][1], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_2), updateSingleHapPanel_->emission_[1][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_2), updateSingleHapPanel_->emission_[1][1], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_3), updateSingleHapPanel_->emission_[2][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_3), updateSingleHapPanel_->emission_[2][1], epsilon3);
    }


    void testUpdateLLK(){
        CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPanel_->updateLLK() );
        CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPlaf_->updateLLK() );
    }


    void testExpectedWsaf(){
        this->updateSingleHapPanel_->segmentStartIndex_ = 0;
        this->updateSingleHapPanel_->nLoci_ = 7;
        this->updateSingleHapPanel_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel_->expectedWsaf0_[0] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.15, this->updateSingleHapPanel_->expectedWsaf0_[1] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.6, this->updateSingleHapPanel_->expectedWsaf0_[2] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel_->expectedWsaf0_[3] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, this->updateSingleHapPanel_->expectedWsaf0_[4] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.55, this->updateSingleHapPanel_->expectedWsaf0_[5] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, this->updateSingleHapPanel_->expectedWsaf0_[6] , epsilon3);
        CPPUNIT_ASSERT_EQUAL ( this->updateSingleHapPanel_->nLoci_, this->updateSingleHapPanel_->expectedWsaf0_.size() );

        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.65, this->updateSingleHapPanel_->expectedWsaf1_[0] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.35, this->updateSingleHapPanel_->expectedWsaf1_[1] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.8, this->updateSingleHapPanel_->expectedWsaf1_[2] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.65, this->updateSingleHapPanel_->expectedWsaf1_[3] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.5, this->updateSingleHapPanel_->expectedWsaf1_[4] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.75, this->updateSingleHapPanel_->expectedWsaf1_[5] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.5, this->updateSingleHapPanel_->expectedWsaf1_[6] , epsilon3);
        CPPUNIT_ASSERT_EQUAL ( this->updateSingleHapPanel_->nLoci_, this->updateSingleHapPanel_->expectedWsaf1_.size() );

        this->updateSingleHapPanel_->strainIndex_ = 3;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.50, this->updateSingleHapPanel_->expectedWsaf0_[0] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.20, this->updateSingleHapPanel_->expectedWsaf0_[1] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel_->expectedWsaf0_[2] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel_->expectedWsaf0_[3] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, this->updateSingleHapPanel_->expectedWsaf0_[4] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.75, this->updateSingleHapPanel_->expectedWsaf0_[5] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, this->updateSingleHapPanel_->expectedWsaf0_[6] , epsilon3);
        CPPUNIT_ASSERT_EQUAL ( this->updateSingleHapPanel_->nLoci_, this->updateSingleHapPanel_->expectedWsaf0_.size() );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.65, this->updateSingleHapPanel_->expectedWsaf1_[0] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.35, this->updateSingleHapPanel_->expectedWsaf1_[1] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.60, this->updateSingleHapPanel_->expectedWsaf1_[2] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.60, this->updateSingleHapPanel_->expectedWsaf1_[3] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel_->expectedWsaf1_[4] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.90, this->updateSingleHapPanel_->expectedWsaf1_[5] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel_->expectedWsaf1_[6] , epsilon3);
        CPPUNIT_ASSERT_EQUAL ( this->updateSingleHapPanel_->nLoci_, this->updateSingleHapPanel_->expectedWsaf1_.size() );

        this->updateSingleHapPlaf_->segmentStartIndex_ = 1;
        this->updateSingleHapPlaf_->nLoci_ = 5;
        this->updateSingleHapPlaf_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPlaf_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.15, this->updateSingleHapPlaf_->expectedWsaf0_[0] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.6, this->updateSingleHapPlaf_->expectedWsaf0_[1] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPlaf_->expectedWsaf0_[2] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, this->updateSingleHapPlaf_->expectedWsaf0_[3] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.55, this->updateSingleHapPlaf_->expectedWsaf0_[4] , epsilon3);
        CPPUNIT_ASSERT_EQUAL ( this->updateSingleHapPlaf_->nLoci_, this->updateSingleHapPlaf_->expectedWsaf0_.size() );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.35, this->updateSingleHapPlaf_->expectedWsaf1_[0] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.8, this->updateSingleHapPlaf_->expectedWsaf1_[1] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.65, this->updateSingleHapPlaf_->expectedWsaf1_[2] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.5, this->updateSingleHapPlaf_->expectedWsaf1_[3] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.75, this->updateSingleHapPlaf_->expectedWsaf1_[4] , epsilon3);
        CPPUNIT_ASSERT_EQUAL ( this->updateSingleHapPlaf_->nLoci_, this->updateSingleHapPlaf_->expectedWsaf1_.size() );
    }


    void testCalcHapLLKs(){
        this->updateSingleHapPanel_->segmentStartIndex_ = 0;
        this->updateSingleHapPanel_->nLoci_ = 7;
        this->updateSingleHapPanel_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel_->calcHapLLKs (this->refCount_, this->altCount_) );

        this->updateSingleHapPlaf_->segmentStartIndex_ = 0;
        this->updateSingleHapPlaf_->nLoci_ = 5;
        this->updateSingleHapPlaf_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPlaf_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPlaf_->calcHapLLKs (this->refCount_, this->altCount_) );
    }


    void testSampleHapIndependently(){
        this->updateSingleHapPlaf_->sampleHapIndependently( this->plaf_ );
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestUpdateSingleHap );
