#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "updateHap.hpp"

class TestUpdatePairHap : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE ( TestUpdatePairHap );
    CPPUNIT_TEST ( testMainConstructor );
    CPPUNIT_TEST ( testComputeMarginalDist );
    CPPUNIT_TEST ( testCalcHapLLKs );
    CPPUNIT_TEST ( testExpectedWsaf );
    CPPUNIT_TEST ( testUpdateLLK );
    CPPUNIT_TEST ( testEmissionProb0 );
    CPPUNIT_TEST_SUITE_END();

  private:
    UpdatePairHap * updatePairHapPanel_;
    UpdatePairHap * updatePairHapPlaf_;
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
    size_t strainIndex1_;
    size_t strainIndex2_;
    bool forbidCopyFromSame_ ;
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
        this->forbidCopyFromSame_ = true;


        this->updatePairHapPlaf_ = new UpdatePairHap( refCount_,
                                          altCount_,
                                          plaf_,
                                          expectedWsaf_,
                                          proportion_,
                                          haplotypes_,
                                          rg_,
                                          segmentStartIndex_,
                                          nLoci_,
                                          NULL,
                                          missCopyProb_, forbidCopyFromSame_,
                                          strainIndex1_, strainIndex2_ );

        this->updatePairHapPanel_ = new UpdatePairHap( refCount_,
                                          altCount_,
                                          plaf_,
                                          expectedWsaf_,
                                          proportion_,
                                          haplotypes_,
                                          rg_,
                                          segmentStartIndex_,
                                          nLoci_,
                                          panel_,
                                          missCopyProb_, forbidCopyFromSame_,
                                          strainIndex1_, strainIndex2_ );
    }


    void tearDown(){
        delete updatePairHapPlaf_;
        delete updatePairHapPanel_;
        delete panel_;
        this->rg_->clearFastFunc();
        delete rg_;
    }


    void testMainConstructor(){ }


    void testComputeMarginalDist(){
        vector < vector < double> > tmpMat;
        tmpMat.push_back( vector <double> ({1.0, 2.0, 3.0}) );
        tmpMat.push_back( vector <double> ({4.0, 5.0, 6.0}) );
        tmpMat.push_back( vector <double> ({7.0, 8.0, 9.0}) );

        vector <double> rowMarg = this->updatePairHapPanel_->computeRowMarginalDist(tmpMat);
        CPPUNIT_ASSERT_EQUAL((size_t)3, rowMarg.size());
        CPPUNIT_ASSERT_EQUAL(6.0,  rowMarg[0]);
        CPPUNIT_ASSERT_EQUAL(15.0, rowMarg[1]);
        CPPUNIT_ASSERT_EQUAL(24.0, rowMarg[2]);

        vector <double> colMarg = this->updatePairHapPanel_->computeColMarginalDist(tmpMat);
        CPPUNIT_ASSERT_EQUAL((size_t)3, colMarg.size());
        CPPUNIT_ASSERT_EQUAL(12.0, colMarg[0]);
        CPPUNIT_ASSERT_EQUAL(15.0, colMarg[1]);
        CPPUNIT_ASSERT_EQUAL(18.0, colMarg[2]);
    }


    void testEmissionProb0 (){
        //double llk0_1 = log(0.05), llk0_2 = log(0.03), llk0_3 = log(0);
        //double llk1_1 = log(0.5), llk1_2 = log(0.0003), llk1_3 = log(1);
        //this->updateSingleHapPanel_->llk0_ = vector <double> ({llk0_1, llk0_2, llk0_3});
        //this->updateSingleHapPanel_->llk1_ = vector <double> ({llk1_1, llk1_2, llk1_3});
        //this->updateSingleHapPanel_->nLoci_ = 3;
        //this->updateSingleHapPanel_->buildEmission ( 0.0 );
        //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05/0.5, updateSingleHapPanel_->emission_[0][0], 0.000000000001);
        //CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel_->emission_[0][1], 0.000000000001);

        //CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel_->emission_[1][0], 0.000000000001);
        //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0003/0.03, updateSingleHapPanel_->emission_[1][1], 0.000000000001);

        //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, updateSingleHapPanel_->emission_[2][0], 0.000000000001);
        //CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel_->emission_[2][1], 0.000000000001);
    }


    void testUpdateLLK(){
        CPPUNIT_ASSERT_NO_THROW ( this->updatePairHapPanel_->updateLLK() );
        CPPUNIT_ASSERT_NO_THROW ( this->updatePairHapPlaf_->updateLLK() );
    }


    void testExpectedWsaf(){
        this->updatePairHapPanel_->segmentStartIndex_ = 0;
        this->updatePairHapPanel_->nLoci_ = 7;
        this->updatePairHapPanel_->strainIndex1_ = 1;
        this->updatePairHapPanel_->strainIndex2_ = 3;
        CPPUNIT_ASSERT_NO_THROW( this->updatePairHapPanel_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
    }


    void testCalcHapLLKs(){
        this->updatePairHapPanel_->segmentStartIndex_ = 0;
        this->updatePairHapPanel_->nLoci_ = 7;
        this->updatePairHapPanel_->strainIndex1_ = 1;
        this->updatePairHapPanel_->strainIndex2_ = 3;
        CPPUNIT_ASSERT_NO_THROW( this->updatePairHapPanel_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updatePairHapPanel_->calcHapLLKs (this->refCount_, this->altCount_) );

        this->updatePairHapPlaf_->segmentStartIndex_ = 0;
        this->updatePairHapPlaf_->nLoci_ = 5;
        this->updatePairHapPlaf_->strainIndex1_ = 1;
        this->updatePairHapPlaf_->strainIndex2_ = 3;
        CPPUNIT_ASSERT_NO_THROW( this->updatePairHapPlaf_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updatePairHapPlaf_->calcHapLLKs (this->refCount_, this->altCount_) );
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestUpdatePairHap );
