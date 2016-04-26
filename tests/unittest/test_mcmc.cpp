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
    CPPUNIT_TEST( testRBernoulli );
    CPPUNIT_TEST( testFindUpdatingStrainSingle );
    CPPUNIT_TEST( testFindUpdatingStrainPair );
    CPPUNIT_TEST( testInitializeHap );
    CPPUNIT_TEST_SUITE_END();

private:
    McmcSample* mcmcSample_;
    PfDeconvIO* pfDeconvIO_;
    Panel* panel_;
    McmcMachinery* mcmcMachinery_;

    void testRBernoulliCore( double prop ){
        size_t nRepeat = 1000000;
        vector <double> rVariables( nRepeat, 0.0 );
        for( vector<double>::iterator it = rVariables.begin(); it != rVariables.end(); ++it) {
            *it = mcmcMachinery_->rBernoulli(prop);
        }
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( prop, sumOfVec(rVariables)/(double)nRepeat, 0.001 );
    }

    void testFindUpdatingStrainSingleCore( size_t kStrain ){
        mcmcMachinery_->kStrain_ = kStrain;
        vector <double> counter ( kStrain, 0.0 );
        size_t nRepeat = 1000000;
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
        size_t nRepeat = 1000000;
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

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestMcmcMachinery );
