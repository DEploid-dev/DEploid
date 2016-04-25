#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "mcmc.hpp"

class TestMcmc: public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestMcmc );

    CPPUNIT_TEST_SUITE_END();

private:


    void setUp() {
    }

    void tearDown() {
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestMcmc );
