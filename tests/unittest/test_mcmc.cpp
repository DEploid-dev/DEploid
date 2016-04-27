#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "mcmc.hpp"

class TestMcmcSample: public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestMcmcSample );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST_SUITE_END();


private:
    McmcSample* mcmcSample_;

public:
    void setUp() {
        mcmcSample_ = new McmcSample();
    }

    void tearDown() {
        delete mcmcSample_;
    }

    void testMainConstructor(){
        McmcSample tmp();
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestMcmcSample );



class TestMcmcMachinery: public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestMcmcMachinery );
    CPPUNIT_TEST( testMainConstructor );
    // Testing for initialization
    CPPUNIT_TEST( testInitializeHap );
    CPPUNIT_TEST( testInitializellk );
    CPPUNIT_TEST( testInitializeProp );
    CPPUNIT_TEST( testInitializeTitre );
    CPPUNIT_TEST( testInitializeExpectedWsaf );
    CPPUNIT_TEST( testCalcMaxIteration );

    CPPUNIT_TEST( testRBernoulli );
    CPPUNIT_TEST( testFindUpdatingStrainSingle );
    CPPUNIT_TEST( testFindUpdatingStrainPair );
    CPPUNIT_TEST( testRunMcmcChain );

    CPPUNIT_TEST_SUITE_END();

private:
    McmcSample* mcmcSample_;
    PfDeconvIO* pfDeconvIO_;
    Panel* panel_;
    McmcMachinery* mcmcMachinery_;
    size_t nRepeat;

    void testRBernoulliCore( double prop ){
        vector <double> rVariables( nRepeat, 0.0 );
        for( vector<double>::iterator it = rVariables.begin(); it != rVariables.end(); ++it) {
            *it = mcmcMachinery_->rBernoulli(prop);
        }
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( prop, sumOfVec(rVariables)/(double)nRepeat, 0.001 );
    }

    void testFindUpdatingStrainSingleCore( size_t kStrain ){
        mcmcMachinery_->kStrain_ = kStrain;
        vector <double> counter ( kStrain, 0.0 );
        for ( size_t i = 0 ; i < nRepeat; i++ ){
            mcmcMachinery_->findUpdatingStrainSingle();
            counter[mcmcMachinery_->strainIndex_] += 1.0 ;
        }
        for ( size_t i = 0 ; i < kStrain; i++){
            CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1.0/(double)kStrain, counter[i]/(double)nRepeat, 0.01 );
        }
    }

    void testFindUpdatingStrainPairCore( size_t kStrain ){
        mcmcMachinery_->kStrain_ = kStrain;
        vector <double> counter ( kStrain, 0.0 );
        for ( size_t i = 0 ; i < nRepeat; i++ ){
            mcmcMachinery_->findUpdatingStrainPair();
            counter[mcmcMachinery_->strainIndex1_] += 1.0 ;
            counter[mcmcMachinery_->strainIndex2_] += 1.0 ;
            CPPUNIT_ASSERT ( mcmcMachinery_->strainIndex1_ != mcmcMachinery_->strainIndex2_ );
        }
        for ( size_t i = 0 ; i < kStrain; i++){
            CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.0/(double)kStrain, counter[i]/(double)nRepeat, 0.01 );
        }
    }


public:
    void setUp() {
        mcmcSample_ = new McmcSample();
        pfDeconvIO_ = new PfDeconvIO();
        panel_ = new Panel();
        mcmcMachinery_ = new McmcMachinery(this->pfDeconvIO_, this->panel_, this->mcmcSample_ );
        nRepeat = 1000000;
    }


    void tearDown() {
        delete mcmcMachinery_;
        delete mcmcSample_;
        delete pfDeconvIO_;
        delete panel_;
    }


    void testMainConstructor(){
        McmcMachinery tmpMcmcMachinery(this->pfDeconvIO_, this->panel_, this->mcmcSample_ );
    }


    void testRBernoulli(){
        this->testRBernoulliCore( 0.23 );
        this->testRBernoulliCore( 0.1 );
        this->testRBernoulliCore( 0.2 );
        this->testRBernoulliCore( 0.4 );
        this->testRBernoulliCore( 0.5 );
        this->testRBernoulliCore( 0.7 );
        this->testRBernoulliCore( 0.8 );
        this->testRBernoulliCore( 0.9 );
    }


    void testFindUpdatingStrainSingle(){
        this->testFindUpdatingStrainSingleCore( 3 );
        this->testFindUpdatingStrainSingleCore( 5 );
        this->testFindUpdatingStrainSingleCore( 6 );
    }


    void testFindUpdatingStrainPair(){
        this->testFindUpdatingStrainPairCore( 3 );
        this->testFindUpdatingStrainPairCore( 5 );
        this->testFindUpdatingStrainPairCore( 6 );
    }


    void testInitializeHap(){
        size_t hapLength = 1000000;
        size_t kStrain = 5;
        size_t constPlaf = 0.3;
        this->mcmcMachinery_->kStrain_ = kStrain;
        this->pfDeconvIO_->plaf_ = vector <double> (hapLength, constPlaf);
        this->mcmcMachinery_->initializeHap();
        for ( size_t i = 0; i < kStrain ; i++ ){
            CPPUNIT_ASSERT_DOUBLES_EQUAL ( constPlaf, sumOfVec(this->mcmcMachinery_->currentHap_[i])/(double)hapLength, 0.001 );
        }
    }


    void testInitializellk(){
        this->mcmcMachinery_->nLoci_ = 1000;
        this->mcmcMachinery_->initializellk();
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.0, sumOfVec(this->mcmcMachinery_->currentLLks_), 0.00001);
        CPPUNIT_ASSERT_EQUAL ( (size_t)1000, this->mcmcMachinery_->currentLLks_.size() );
    }


    void testInitializeProp(){
        // TODO
    }


    void testInitializeTitre(){
        // TODO
    }


    void testInitializeExpectedWsaf(){
        // TODO
    }


    void testDeltaLLKs(){
        this->mcmcMachinery_->currentLLks_ = vector < double > ( { 0.1, 0.3, 0.2, 0.4 } );
        vector <double> newllks ({0.3,0.2, 0.5, 0.6});
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.6, this->mcmcMachinery_->deltaLLKs( newllks ) , 0.00001);
    }


    void testCalcMaxIteration(){
        CPPUNIT_ASSERT_EQUAL( 0.5, this->mcmcMachinery_->burnIn_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)5, this->mcmcMachinery_->McmcMachineryRate_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)8001, this->mcmcMachinery_->maxIteration_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)4000, this->mcmcMachinery_->mcmcThresh_ );

        this->mcmcMachinery_->calcMaxIteration(1000, 5, 0.5);
        CPPUNIT_ASSERT_EQUAL( 0.5, this->mcmcMachinery_->burnIn_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)5, this->mcmcMachinery_->McmcMachineryRate_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)10001, this->mcmcMachinery_->maxIteration_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)5000, this->mcmcMachinery_->mcmcThresh_ );

        this->mcmcMachinery_->calcMaxIteration(10, 3, 0.3);
        CPPUNIT_ASSERT_EQUAL( 0.3, this->mcmcMachinery_->burnIn_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)3, this->mcmcMachinery_->McmcMachineryRate_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)44, this->mcmcMachinery_->maxIteration_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)13, this->mcmcMachinery_->mcmcThresh_ );

        this->mcmcMachinery_->calcMaxIteration(10, 7, 0.4);
        CPPUNIT_ASSERT_EQUAL( 0.4, this->mcmcMachinery_->burnIn_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)7, this->mcmcMachinery_->McmcMachineryRate_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)118, this->mcmcMachinery_->maxIteration_ );
        CPPUNIT_ASSERT_EQUAL( (size_t)47, this->mcmcMachinery_->mcmcThresh_ );
    }


    void testRunMcmcChain (){
        this->mcmcMachinery_->calcMaxIteration(nRepeat, 1, 0.5);
        this->mcmcMachinery_->runMcmcChain();

        CPPUNIT_ASSERT_EQUAL( nRepeat, this->mcmcMachinery_->mcmcSample_->proportion.size() );
        CPPUNIT_ASSERT_EQUAL( nRepeat, this->mcmcMachinery_->mcmcSample_->sumLLKs.size() );
        CPPUNIT_ASSERT_EQUAL( nRepeat, this->mcmcMachinery_->mcmcSample_->moves.size() );
        vector <size_t> counter (3, 0);
        for ( size_t i = 0; i < nRepeat; i++){
            counter[this->mcmcMachinery_->mcmcSample_->moves[i]] += 1;
        }
        for ( size_t i = 0; i < 3; i++ ){
            CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.333333, (double)counter[i]/(double)nRepeat , 0.001);
        }
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestMcmcMachinery );
