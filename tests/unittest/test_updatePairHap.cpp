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
    CPPUNIT_TEST ( testSampleHapIndependently );
    CPPUNIT_TEST ( testCore );

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
    double epsilon2;
    size_t nRepeat;

  public:
    void setUp(){
        nRepeat = 1000000;
        epsilon2 = 0.001;
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
        double llk00_1 = log(0.05), llk00_2 = log(0.03),   llk00_3 = log(0);
        double llk01_1 = log(0.5),  llk01_2 = log(0.0003), llk01_3 = log(1);
        double llk10_1 = log(0.35), llk10_2 = log(0.3),    llk10_3 = log(0.3);
        double llk11_1 = log(0.65), llk11_2 = log(0.43),   llk11_3 = log(0.1);
        this->updatePairHapPanel_->llk00_ = vector <double> ({llk00_1, llk00_2, llk00_3});
        this->updatePairHapPanel_->llk01_ = vector <double> ({llk01_1, llk01_2, llk01_3});
        this->updatePairHapPanel_->llk10_ = vector <double> ({llk10_1, llk10_2, llk10_3});
        this->updatePairHapPanel_->llk11_ = vector <double> ({llk11_1, llk11_2, llk11_3});
        this->updatePairHapPanel_->nLoci_ = 3;
        CPPUNIT_ASSERT_NO_THROW ( this->updatePairHapPanel_->buildEmission ( 0.0 ) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk00_1-llk11_1), updatePairHapPanel_->emission_[0][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk01_1-llk11_1), updatePairHapPanel_->emission_[0][1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk10_1-llk11_1), updatePairHapPanel_->emission_[0][2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk11_1-llk11_1), updatePairHapPanel_->emission_[0][3], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk00_2-llk11_2), updatePairHapPanel_->emission_[1][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk01_2-llk11_2), updatePairHapPanel_->emission_[1][1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk10_2-llk11_2), updatePairHapPanel_->emission_[1][2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk11_2-llk11_2), updatePairHapPanel_->emission_[1][3], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk00_3-llk01_3), updatePairHapPanel_->emission_[2][0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk01_3-llk01_3), updatePairHapPanel_->emission_[2][1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk10_3-llk01_3), updatePairHapPanel_->emission_[2][2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(llk11_3-llk01_3), updatePairHapPanel_->emission_[2][3], epsilon3);

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


    void testSampleHapIndependently(){
        this->updatePairHapPanel_->segmentStartIndex_ = 0;
        this->updatePairHapPanel_->nLoci_ = 7;
        this->updatePairHapPanel_->strainIndex1_ = 1;
        this->updatePairHapPanel_->strainIndex2_ = 3;
        CPPUNIT_ASSERT_NO_THROW( this->updatePairHapPanel_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updatePairHapPanel_->calcHapLLKs (this->refCount_, this->altCount_) );


        //for ( size_t i = 0 ; i < 7;i++){
                //double tmpMax = max_value ( vector <double> ( {this->updatePairHapPanel_->llk00_[i],
                                                               //this->updatePairHapPanel_->llk01_[i],
                                                               //this->updatePairHapPanel_->llk10_[i],
                                                               //this->updatePairHapPanel_->llk11_[i]} ) );
                //cout << exp(this->updatePairHapPanel_->llk00_[i]-tmpMax) <<" "
                     //<< exp(this->updatePairHapPanel_->llk01_[i]-tmpMax) <<" "
                     //<< exp(this->updatePairHapPanel_->llk10_[i]-tmpMax) <<" "
                     //<< exp(this->updatePairHapPanel_->llk11_[i]-tmpMax) << endl;
        //}

//1 1.50559e-06 9.59254e-09 2.58807e-16
//5.17896e-31 1.46301e-13 9.15637e-10 1
//1 0.476256 0.134298 6.73633e-05
//1 0.572198 0.219935 0.000642945
//1 4.49961e-06 3.83615e-08 2.1567e-15
//3.58261e-07 0.00196058 0.018196 1
//0.992688 1 0.412035 0.00204737

        this->plaf_ = vector < double > (this->updatePairHapPanel_->nLoci_, 0.5);
        vector < vector <int> > counter;
        for ( size_t i = 0; i < this->updatePairHapPanel_->nLoci_; i++ ){
            counter.push_back ( vector <int> (4, 0) );
        }

        for ( size_t i = 0; i < nRepeat; i++ ){
            this->updatePairHapPanel_->hap1_.clear();
            this->updatePairHapPanel_->hap2_.clear();
            CPPUNIT_ASSERT_NO_THROW ( this->updatePairHapPanel_->sampleHapIndependently( this->plaf_ ) );
            for ( size_t j = 0; j < this->updatePairHapPanel_->nLoci_; j++ ){
                counter[j][(size_t)(this->updatePairHapPanel_->hap1_[j]*2+this->updatePairHapPanel_->hap2_[j])] ++ ;
            }
        }

        CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, double(counter[0][0])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[0][1])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[0][2])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[0][3])/(double)nRepeat, epsilon2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[1][0])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[1][1])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[1][2])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, double(counter[1][3])/(double)nRepeat, epsilon2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.621, double(counter[2][0])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.296, double(counter[2][1])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.083, double(counter[2][2])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[2][3])/(double)nRepeat, epsilon2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.558, double(counter[3][0])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.319, double(counter[3][1])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.123, double(counter[3][2])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[3][3])/(double)nRepeat, epsilon2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, double(counter[4][0])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[4][1])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[4][2])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[4][3])/(double)nRepeat, epsilon2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[5][0])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.002, double(counter[5][1])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.018, double(counter[5][2])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.98, double(counter[5][3])/(double)nRepeat, epsilon2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.412, double(counter[6][0])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.415, double(counter[6][1])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.171, double(counter[6][2])/(double)nRepeat, epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL (0.001, double(counter[6][3])/(double)nRepeat, epsilon2);


        this->plaf_ = vector < double > (this->updatePairHapPanel_->nLoci_, 0.0);
        counter.clear();
        for ( size_t i = 0; i < this->updatePairHapPanel_->nLoci_; i++ ){
            counter.push_back ( vector <int> (4, 0) );
        }

        for ( size_t i = 0; i < nRepeat; i++ ){
            this->updatePairHapPanel_->hap1_.clear();
            this->updatePairHapPanel_->hap2_.clear();
            CPPUNIT_ASSERT_NO_THROW ( this->updatePairHapPanel_->sampleHapIndependently( this->plaf_ ) );
            for ( size_t j = 0; j < this->updatePairHapPanel_->nLoci_; j++ ){
                counter[j][(size_t)(this->updatePairHapPanel_->hap1_[j]*2+this->updatePairHapPanel_->hap2_[j])] ++ ;
            }
        }
        for ( size_t j = 0; j < this->updatePairHapPanel_->nLoci_; j++ ){
            CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, double(counter[j][0])/(double)nRepeat, epsilon2);
            CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[j][1])/(double)nRepeat, epsilon2);
            CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[j][2])/(double)nRepeat, epsilon2);
            CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[j][3])/(double)nRepeat, epsilon2);
        }

        this->plaf_ = vector < double > (this->updatePairHapPanel_->nLoci_, 1.0);
        counter.clear();
        for ( size_t i = 0; i < this->updatePairHapPanel_->nLoci_; i++ ){
            counter.push_back ( vector <int> (4, 0) );
        }

        for ( size_t i = 0; i < nRepeat; i++ ){
            this->updatePairHapPanel_->hap1_.clear();
            this->updatePairHapPanel_->hap2_.clear();
            CPPUNIT_ASSERT_NO_THROW ( this->updatePairHapPanel_->sampleHapIndependently( this->plaf_ ) );
            for ( size_t j = 0; j < this->updatePairHapPanel_->nLoci_; j++ ){
                counter[j][(size_t)(this->updatePairHapPanel_->hap1_[j]*2+this->updatePairHapPanel_->hap2_[j])] ++ ;
            }
        }
        for ( size_t j = 0; j < this->updatePairHapPanel_->nLoci_; j++ ){
            CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[j][0])/(double)nRepeat, epsilon2);
            CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[j][1])/(double)nRepeat, epsilon2);
            CPPUNIT_ASSERT_DOUBLES_EQUAL (0.0, double(counter[j][2])/(double)nRepeat, epsilon2);
            CPPUNIT_ASSERT_DOUBLES_EQUAL (1.0, double(counter[j][3])/(double)nRepeat, epsilon2);
        }
    }


    void testUpdateLLK(){
        CPPUNIT_ASSERT_NO_THROW ( this->updatePairHapPanel_->updateLLK() );

        this->updatePairHapPlaf_->segmentStartIndex_ = 0;
        this->updatePairHapPlaf_->nLoci_ = 7;
        this->updatePairHapPanel_->strainIndex1_ = 1;
        this->updatePairHapPanel_->strainIndex2_ = 3;
        CPPUNIT_ASSERT_NO_THROW( this->updatePairHapPlaf_->calcExpectedWsaf( this->expectedWsaf_, this->proportion_, this->haplotypes_) );
        CPPUNIT_ASSERT_NO_THROW( this->updatePairHapPlaf_->calcHapLLKs (this->refCount_, this->altCount_) );
        this->plaf_ = vector < double > (this->updatePairHapPlaf_->nLoci_, 0.5);
        CPPUNIT_ASSERT_NO_THROW ( this->updatePairHapPlaf_->sampleHapIndependently( this->plaf_ ) );
        CPPUNIT_ASSERT_NO_THROW ( this->updatePairHapPlaf_->updateLLK() );
    }


    void testCore(){
        this->plaf_ = vector < double > (this->updatePairHapPlaf_->nLoci_, 0.5);
        CPPUNIT_ASSERT_NO_THROW( this->updatePairHapPlaf_->core( this->refCount_,
                                                                 this->altCount_,
                                                                 this->plaf_,
                                                                 this->expectedWsaf_,
                                                                 this->proportion_,
                                                                 this->haplotypes_ ) );
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestUpdatePairHap );
