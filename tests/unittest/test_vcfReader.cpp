#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "src/vcfReader.hpp"

class TestVCF : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestVCF );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST_SUITE_END();

  private:
    VcfReader* vcf_;
    VcfReader* vcfGz_;

  public:
    void setUp() {
        this->vcf_ = new VcfReader ( "data/testData/PG0390-C.test.vcf" );
        this->vcfGz_ = new VcfReader ( "data/testData/PG0390-C.test.vcf.gz" );
    }

    void tearDown(){
        delete this->vcf_;
        delete this->vcfGz_;
    }

    void testMainConstructor(){
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestVCF );
