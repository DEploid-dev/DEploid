#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "vcfReader.hpp"

class TestVCF : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestVCF );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST_SUITE_END();

  private:
    VcfReader* vcf_;

  public:
    void setUp() {
        this->vcf_ = new VcfReader ( "data/testData/PG0390-C.test.vcf" );
    }

    void tearDown(){
        delete this->vcf_;
    }

    void testMainConstructor(){
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestVCF );
