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
        delete rg_;
    }


    void testMainConstructor(){ }


    void testVirtualFunctions(){
        CPPUNIT_ASSERT_THROW ( this->updateHapPlaf_->core ( refCount_, altCount_, plaf_, expectedWsaf_, proportion_, haplotypes_), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPlaf_->calcExpectedWsaf(expectedWsaf_, proportion_, haplotypes_ ), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPlaf_->calcHapLLKs(refCount_, altCount_), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPlaf_->buildEmission(0.1), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPlaf_->samplePaths(), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPlaf_->addMissCopying(0.1), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPlaf_->updateLLK(), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPlaf_->sampleHapIndependently(plaf_), VirtualFunctionShouldNotBeCalled );

        CPPUNIT_ASSERT_THROW ( this->updateHapPanel_->core ( refCount_, altCount_, plaf_, expectedWsaf_, proportion_, haplotypes_ ), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPanel_->calcExpectedWsaf(expectedWsaf_, proportion_, haplotypes_ ), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPanel_->calcHapLLKs(refCount_, altCount_), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPanel_->buildEmission(0.1), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPanel_->samplePaths(), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPanel_->addMissCopying(0.1), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPanel_->updateLLK(), VirtualFunctionShouldNotBeCalled );
        CPPUNIT_ASSERT_THROW ( this->updateHapPanel_->sampleHapIndependently(plaf_), VirtualFunctionShouldNotBeCalled );
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestUpdateHap );


class TestUpdateSingleHap : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestUpdateSingleHap );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST(testFwdBwd);
    //CPPUNIT_TEST ( testEmissionProb0 );
    //CPPUNIT_TEST ( testEmissionBasicVersion0 );
    //CPPUNIT_TEST ( testUpdateLLK );
    //CPPUNIT_TEST ( testExpectedWsaf );
    //CPPUNIT_TEST ( testCalcHapLLKs );
    //CPPUNIT_TEST ( testSampleHapIndependently );
    //CPPUNIT_TEST ( testCore );
    //CPPUNIT_TEST ( testCalcFwdProbsNoRecomb );
    CPPUNIT_TEST_SUITE_END();

  private:
    UpdateSingleHap * updateSingleHapPanel1_;
    UpdateSingleHap * updateSingleHapPanel2_;
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
    Panel* panel1_;
    Panel* panel2_;
    double missCopyProb_;
    size_t strainIndex_;
    double epsilon1;
    double epsilon2;
    double epsilon3;
    size_t nRepeat;


  public:

    void setUp(){
        nRepeat = 1000000;
        epsilon1 = 0.001;
        epsilon2 = 0.000001;
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
        this->panel1_ = new Panel();
        this->panel1_->buildExamplePanel1();

        this->panel2_ = new Panel();
        this->panel2_->buildExamplePanel2();

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

        this->updateSingleHapPanel1_ = new UpdateSingleHap( refCount_,
                                          altCount_,
                                          plaf_,
                                          expectedWsaf_,
                                          proportion_,
                                          haplotypes_,
                                          rg_,
                                          segmentStartIndex_,
                                          nLoci_,
                                          panel1_,
                                          missCopyProb_,
                                          strainIndex_ );

        this->updateSingleHapPanel2_ = new UpdateSingleHap( refCount_,
                                          altCount_,
                                          plaf_,
                                          expectedWsaf_,
                                          proportion_,
                                          haplotypes_,
                                          rg_,
                                          segmentStartIndex_,
                                          nLoci_,
                                          panel2_,
                                          missCopyProb_,
                                          strainIndex_ );
}


    void tearDown(){
        delete updateSingleHapPlaf_;
        delete updateSingleHapPanel1_;
        delete updateSingleHapPanel2_;
        delete panel1_;
        delete panel2_;
        delete rg_;
    }


    void testMainConstructor(){ }


void testFwdBwd(){
        CPPUNIT_ASSERT_NO_THROW(this->panel1_->computeRecombProbs( 15000.0, 10.0, true, 0, false));

        this->updateSingleHapPanel1_->segmentStartIndex_ = 0;
        this->updateSingleHapPanel1_->nLoci_ = 7;
        this->updateSingleHapPanel1_->strainIndex_ = 0;

        //CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        //CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->calcHapLLKs (this->refCount_, this->altCount_) );

        this->updateSingleHapPanel1_->llk0_ = vector <double> (7, log(0.1));
        this->updateSingleHapPanel1_->llk1_ = vector <double> (7, log(0.9));

        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->buildEmission ( 0.0 ) );
        //CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->buildEmission ( 0.01 ) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->calcFwdBwdProbs() );

cout <<endl;
        for ( vector<double> pRow : this->updateSingleHapPanel1_->emission_){
            for (double p : pRow){
                cout<<p<< " ";
            }
            cout <<endl;
        }



cout <<endl;
        for ( vector<double> pRow : this->updateSingleHapPanel1_->fwdProbs_){
            for (double p : pRow){
                cout<<p<< " ";
            }
            cout <<endl;
        }
cout <<endl;
        for ( vector<double> pRow : this->updateSingleHapPanel1_->bwdProbs_){
            for (double p : pRow){
                cout<<p<< " ";
            }
            cout <<endl;
        }
cout <<endl;
        for ( vector<double> pRow : this->updateSingleHapPanel1_->fwdBwdProbs_){
            for (double p : pRow){
                cout<<p<< " ";
            }
            cout <<endl;
        }

        //cout <<this->updateSingleHapPanel1_->fwdProbs_.size() << endl;

}

    void testEmissionProb0 (){
        double llk0_1 = log(0.05), llk0_2 = log(0.03), llk0_3 = log(0);
        double llk1_1 = log(0.5), llk1_2 = log(0.0003), llk1_3 = log(1);
        this->updateSingleHapPanel1_->llk0_ = vector <double> ({llk0_1, llk0_2, llk0_3});
        this->updateSingleHapPanel1_->llk1_ = vector <double> ({llk1_1, llk1_2, llk1_3});
        this->updateSingleHapPanel1_->nLoci_ = 3;
        CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPanel1_->buildEmission ( 0.0 ) );
        CPPUNIT_ASSERT_EQUAL ( (size_t)2, this->updateSingleHapPanel1_->emission_[0].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)2, this->updateSingleHapPanel1_->emission_[1].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)2, this->updateSingleHapPanel1_->emission_[2].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)2, this->updateSingleHapPanel1_->emission_.back().size() );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_1-llk1_1), updateSingleHapPanel1_->emission_[0][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05/0.5, updateSingleHapPanel1_->emission_[0][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel1_->emission_[0][1], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel1_->emission_[1][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0003/0.03, updateSingleHapPanel1_->emission_[1][1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_2-llk0_2), updateSingleHapPanel1_->emission_[1][1], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_3-llk1_3), updateSingleHapPanel1_->emission_[2][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, updateSingleHapPanel1_->emission_[2][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, updateSingleHapPanel1_->emission_[2][1], epsilon3);
    }


    void testEmissionBasicVersion0 (){
        double llk0_1 = log(0.05), llk0_2 = log(0.03), llk0_3 = log(0);
        double llk1_1 = log(0.5), llk1_2 = log(0.0003), llk1_3 = log(1);
        this->updateSingleHapPanel1_->llk0_ = vector <double> ({llk0_1, llk0_2, llk0_3});
        this->updateSingleHapPanel1_->llk1_ = vector <double> ({llk1_1, llk1_2, llk1_3});
        this->updateSingleHapPanel1_->nLoci_ = 3;
        CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPanel1_->buildEmissionBasicVersion ( 0.0 ) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_1), updateSingleHapPanel1_->emission_[0][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_1), updateSingleHapPanel1_->emission_[0][1], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_2), updateSingleHapPanel1_->emission_[1][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_2), updateSingleHapPanel1_->emission_[1][1], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk0_3), updateSingleHapPanel1_->emission_[2][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk1_3), updateSingleHapPanel1_->emission_[2][1], epsilon3);
    }


    void testExpectedWsaf(){
        this->updateSingleHapPanel1_->segmentStartIndex_ = 0;
        this->updateSingleHapPanel1_->nLoci_ = 7;
        this->updateSingleHapPanel1_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel1_->expectedWsaf0_[0] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.15, this->updateSingleHapPanel1_->expectedWsaf0_[1] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.6, this->updateSingleHapPanel1_->expectedWsaf0_[2] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel1_->expectedWsaf0_[3] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, this->updateSingleHapPanel1_->expectedWsaf0_[4] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.55, this->updateSingleHapPanel1_->expectedWsaf0_[5] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, this->updateSingleHapPanel1_->expectedWsaf0_[6] , epsilon3);
        CPPUNIT_ASSERT_EQUAL ( this->updateSingleHapPanel1_->nLoci_, this->updateSingleHapPanel1_->expectedWsaf0_.size() );

        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.65, this->updateSingleHapPanel1_->expectedWsaf1_[0] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.35, this->updateSingleHapPanel1_->expectedWsaf1_[1] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.8, this->updateSingleHapPanel1_->expectedWsaf1_[2] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.65, this->updateSingleHapPanel1_->expectedWsaf1_[3] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.5, this->updateSingleHapPanel1_->expectedWsaf1_[4] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.75, this->updateSingleHapPanel1_->expectedWsaf1_[5] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.5, this->updateSingleHapPanel1_->expectedWsaf1_[6] , epsilon3);
        CPPUNIT_ASSERT_EQUAL ( this->updateSingleHapPanel1_->nLoci_, this->updateSingleHapPanel1_->expectedWsaf1_.size() );

        this->updateSingleHapPanel1_->strainIndex_ = 3;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.50, this->updateSingleHapPanel1_->expectedWsaf0_[0] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.20, this->updateSingleHapPanel1_->expectedWsaf0_[1] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel1_->expectedWsaf0_[2] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel1_->expectedWsaf0_[3] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, this->updateSingleHapPanel1_->expectedWsaf0_[4] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.75, this->updateSingleHapPanel1_->expectedWsaf0_[5] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.3, this->updateSingleHapPanel1_->expectedWsaf0_[6] , epsilon3);
        CPPUNIT_ASSERT_EQUAL ( this->updateSingleHapPanel1_->nLoci_, this->updateSingleHapPanel1_->expectedWsaf0_.size() );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.65, this->updateSingleHapPanel1_->expectedWsaf1_[0] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.35, this->updateSingleHapPanel1_->expectedWsaf1_[1] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.60, this->updateSingleHapPanel1_->expectedWsaf1_[2] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.60, this->updateSingleHapPanel1_->expectedWsaf1_[3] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel1_->expectedWsaf1_[4] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.90, this->updateSingleHapPanel1_->expectedWsaf1_[5] , epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.45, this->updateSingleHapPanel1_->expectedWsaf1_[6] , epsilon3);
        CPPUNIT_ASSERT_EQUAL ( this->updateSingleHapPanel1_->nLoci_, this->updateSingleHapPanel1_->expectedWsaf1_.size() );

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
        this->updateSingleHapPanel1_->segmentStartIndex_ = 0;
        this->updateSingleHapPanel1_->nLoci_ = 7;
        this->updateSingleHapPanel1_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->calcHapLLKs (this->refCount_, this->altCount_) );

        this->updateSingleHapPlaf_->segmentStartIndex_ = 0;
        this->updateSingleHapPlaf_->nLoci_ = 5;
        this->updateSingleHapPlaf_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPlaf_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPlaf_->calcHapLLKs (this->refCount_, this->altCount_) );
    }


    void testSampleHapIndependently(){
        this->updateSingleHapPlaf_->segmentStartIndex_ = 0;
        this->updateSingleHapPlaf_->nLoci_ = 7;
        this->updateSingleHapPlaf_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPlaf_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPlaf_->calcHapLLKs (this->refCount_, this->altCount_) );
        //for ( size_t i = 0 ; i < 7;i++){
                //double tmpMax = max(this->updateSingleHapPlaf_->llk0_[i], this->updateSingleHapPlaf_->llk1_[i]);
                //cout << exp(this->updateSingleHapPlaf_->llk0_[i]-tmpMax) <<" " << exp(this->updateSingleHapPlaf_->llk1_[i]-tmpMax) << endl;
        //}
//1 1.71897e-10
//1.46301e-13 1
//1 0.000141443
//1 0.219935
//1 3.83615e-08
//1.9689e-05 1
//1 0.41507

        this->plaf_ = vector < double > (this->updateSingleHapPlaf_->nLoci_, 0.5);
        vector < vector <int> > counter;
        for ( size_t i = 0; i < this->updateSingleHapPlaf_->nLoci_; i++ ){
            counter.push_back ( vector <int> (2, 0) );
        }

        for ( size_t i = 0; i < nRepeat; i++ ){
            this->updateSingleHapPlaf_->hap_.clear();
            CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPlaf_->sampleHapIndependently( this->plaf_ ) );
            for ( size_t j = 0; j < this->updateSingleHapPlaf_->nLoci_; j++ ){
                counter[j][(size_t)this->updateSingleHapPlaf_->hap_[j]] ++ ;
            }
        }

        CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, double(counter[0][0])/(double)nRepeat, epsilon1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[0][1])/(double)nRepeat, epsilon1);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[1][0])/(double)nRepeat, epsilon1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, double(counter[1][1])/(double)nRepeat, epsilon1);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.999, double(counter[2][0])/(double)nRepeat, epsilon1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.001, double(counter[2][1])/(double)nRepeat, epsilon1);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.82, double(counter[3][0])/(double)nRepeat, epsilon1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.18, double(counter[3][1])/(double)nRepeat, epsilon1);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, double(counter[4][0])/(double)nRepeat, epsilon1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[4][1])/(double)nRepeat, epsilon1);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[5][0])/(double)nRepeat, epsilon1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, double(counter[5][1])/(double)nRepeat, epsilon1);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.707, double(counter[6][0])/(double)nRepeat, epsilon1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.293, double(counter[6][1])/(double)nRepeat, epsilon1);

        this->plaf_ = vector < double > (this->updateSingleHapPlaf_->nLoci_, 0.0);
        counter.clear();
        for ( size_t i = 0; i < this->updateSingleHapPlaf_->nLoci_; i++ ){
            counter.push_back ( vector <int> (2, 0) );
        }

        for ( size_t i = 0; i < nRepeat; i++ ){
            this->updateSingleHapPlaf_->hap_.clear();
            CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPlaf_->sampleHapIndependently( this->plaf_ ) );
            for ( size_t j = 0; j < this->updateSingleHapPlaf_->nLoci_; j++ ){
                counter[j][(size_t)this->updateSingleHapPlaf_->hap_[j]] ++ ;
            }
        }
        for ( size_t j = 0; j < this->updateSingleHapPlaf_->nLoci_; j++ ){
            CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, double(counter[j][0])/(double)nRepeat, epsilon1);
            CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[j][1])/(double)nRepeat, epsilon1);
        }


        this->plaf_ = vector < double > (this->updateSingleHapPlaf_->nLoci_, 1.0);
        counter.clear();
        for ( size_t i = 0; i < this->updateSingleHapPlaf_->nLoci_; i++ ){
            counter.push_back ( vector <int> (2, 0) );
        }

        for ( size_t i = 0; i < nRepeat; i++ ){
            this->updateSingleHapPlaf_->hap_.clear();
            CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPlaf_->sampleHapIndependently( this->plaf_ ) );
            for ( size_t j = 0; j < this->updateSingleHapPlaf_->nLoci_; j++ ){
                counter[j][(size_t)this->updateSingleHapPlaf_->hap_[j]] ++ ;
            }
        }
        for ( size_t j = 0; j < this->updateSingleHapPlaf_->nLoci_; j++ ){
            CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[j][0])/(double)nRepeat, epsilon1);
            CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, double(counter[j][1])/(double)nRepeat, epsilon1);
        }

    }


    void testUpdateLLK(){
        CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPanel1_->updateLLK() );

        this->updateSingleHapPlaf_->segmentStartIndex_ = 0;
        this->updateSingleHapPlaf_->nLoci_ = 7;
        this->updateSingleHapPlaf_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPlaf_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPlaf_->calcHapLLKs (this->refCount_, this->altCount_) );
        this->plaf_ = vector < double > (this->updateSingleHapPlaf_->nLoci_, 0.5);
        CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPlaf_->sampleHapIndependently( this->plaf_ ) );
        CPPUNIT_ASSERT_NO_THROW ( this->updateSingleHapPlaf_->updateLLK() );
    }


    void testCore(){
        this->plaf_ = vector < double > (this->updateSingleHapPlaf_->nLoci_, 0.5);
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPlaf_->core( this->refCount_,
                                                                   this->altCount_,
                                                                   this->plaf_,
                                                                   this->expectedWsaf_,
                                                                   this->proportion_,
                                                                   this->haplotypes_ ) );

        CPPUNIT_ASSERT_NO_THROW(this->panel1_->computeRecombProbs( 15000.0, 10.0, true, 0, false));
        this->updateSingleHapPanel1_->segmentStartIndex_ = 0;
        this->updateSingleHapPanel1_->nLoci_ = 7;
        this->updateSingleHapPanel1_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->core( this->refCount_,
                                                                   this->altCount_,
                                                                   this->plaf_,
                                                                   this->expectedWsaf_,
                                                                   this->proportion_,
                                                                   this->haplotypes_ ) );

        CPPUNIT_ASSERT_NO_THROW(this->panel2_->computeRecombProbs( 15000.0, 10.0, true, 0, false));
        this->updateSingleHapPanel2_->segmentStartIndex_ = 1;
        this->updateSingleHapPanel2_->nLoci_ = 5;
        this->updateSingleHapPanel2_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel2_->core( this->refCount_,
                                                                   this->altCount_,
                                                                   this->plaf_,
                                                                   this->expectedWsaf_,
                                                                   this->proportion_,
                                                                   this->haplotypes_ ) );
    }


    void testCalcFwdProbsNoRecomb(){
        CPPUNIT_ASSERT_NO_THROW(this->panel1_->computeRecombProbs( 15000.0, 10.0, true, 0, false));

        this->updateSingleHapPanel1_->segmentStartIndex_ = 0;
        this->updateSingleHapPanel1_->nLoci_ = 7;
        this->updateSingleHapPanel1_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->calcHapLLKs (this->refCount_, this->altCount_) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->buildEmission ( 0.0 ) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel1_->calcFwdProbs() );
//StrainHap: 1,1,0,0,0,1,0
//this->refCount_ = vector <double> (  { 100, 10 , 50 , 30 , 100, 7, 50 } );
//this->altCount_ = vector <double> (  { 2,  100 , 50 , 30 , 5,  70, 30 } );
    //this->content_.push_back( vector <double> ({0,0,0,1}) );
    //this->content_.push_back( vector <double> ({0,0,0,1}) );
    //this->content_.push_back( vector <double> ({0,0,0,1}) );
    //this->content_.push_back( vector <double> ({0,0,0,1}) );
    //this->content_.push_back( vector <double> ({0,1,1,0}) );
    //this->content_.push_back( vector <double> ({0,0,1,0}) );
    //this->content_.push_back( vector <double> ({0,0,1,0}) );
        // At the first site, Ref count 100, greater larger than alt, emmision to 0 with prob 1, to 1 with prob 0, thus, first three panel hap are equally likely
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.3333333, this->updateSingleHapPanel1_->fwdProbs_[0][0], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, this->updateSingleHapPanel1_->fwdProbs_[0][3], epsilon2);

        // At the second site, Ref count 10, far less than alt, emmision to 0 with prob 0, to 1 with prob 1
        // loci at postion [0]-[3], first three haps should all have equal probabilities.
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel1_->fwdProbs_[0][0], this->updateSingleHapPanel1_->fwdProbs_[0][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel1_->fwdProbs_[0][2], this->updateSingleHapPanel1_->fwdProbs_[0][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel1_->fwdProbs_[1][0], this->updateSingleHapPanel1_->fwdProbs_[1][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel1_->fwdProbs_[1][2], this->updateSingleHapPanel1_->fwdProbs_[1][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel1_->fwdProbs_[2][0], this->updateSingleHapPanel1_->fwdProbs_[2][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel1_->fwdProbs_[2][2], this->updateSingleHapPanel1_->fwdProbs_[2][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel1_->fwdProbs_[3][0], this->updateSingleHapPanel1_->fwdProbs_[3][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel1_->fwdProbs_[3][2], this->updateSingleHapPanel1_->fwdProbs_[3][1], epsilon2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel1_->fwdProbs_[4][2], this->updateSingleHapPanel1_->fwdProbs_[4][1], epsilon2);

        CPPUNIT_ASSERT_NO_THROW(this->panel2_->computeRecombProbs( 15000.0, 10.0, true, 0, false));

        this->updateSingleHapPanel2_->segmentStartIndex_ = 1;
        this->updateSingleHapPanel2_->nLoci_ = 5;
        this->updateSingleHapPanel2_->strainIndex_ = 1;
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel2_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel2_->calcHapLLKs (this->refCount_, this->altCount_) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel2_->buildEmission ( 0.0 ) );
        CPPUNIT_ASSERT_NO_THROW( this->updateSingleHapPanel2_->calcFwdProbs() );
    //this->content_.push_back( vector <double> ({0,0,0,1}) );
    //this->content_.push_back( vector <double> ({0,0,0,1}) );
    //this->content_.push_back( vector <double> ({0,0,0,1}) );
    //this->content_.push_back( vector <double> ({0,1,1,0}) );
    //this->content_.push_back( vector <double> ({0,0,1,0}) );

        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel2_->fwdProbs_[0][0], this->updateSingleHapPanel2_->fwdProbs_[0][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel2_->fwdProbs_[0][2], this->updateSingleHapPanel2_->fwdProbs_[0][1], epsilon2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, this->updateSingleHapPanel2_->fwdProbs_[0][3], epsilon2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel2_->fwdProbs_[1][0], this->updateSingleHapPanel2_->fwdProbs_[1][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel2_->fwdProbs_[1][2], this->updateSingleHapPanel2_->fwdProbs_[1][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel2_->fwdProbs_[2][0], this->updateSingleHapPanel2_->fwdProbs_[2][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel2_->fwdProbs_[2][2], this->updateSingleHapPanel2_->fwdProbs_[2][1], epsilon2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (this->updateSingleHapPanel2_->fwdProbs_[3][2], this->updateSingleHapPanel2_->fwdProbs_[3][1], epsilon2);

    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestUpdateSingleHap );
