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

    CPPUNIT_TEST_SUITE_END();


private:
    McmcSample* mcmcSample_;
    PfDeconvIO* pfDeconvIO_;
    Panel* panel_;

public:
    void setUp() {
        mcmcSample_ = new McmcSample();
        pfDeconvIO_ = new PfDeconvIO();
        panel_ = new Panel();
    }

    void tearDown() {
        delete mcmcSample_;
        delete pfDeconvIO_;
        delete panel_;
    }

    void testMainConstructor(){
        McmcMachinery tmpMcmcMachinery(this->pfDeconvIO_, this->panel_, this->mcmcSample_ );
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestMcmcMachinery );
