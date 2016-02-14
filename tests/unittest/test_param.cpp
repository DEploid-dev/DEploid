#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../../src/param.hpp"

class TestInput : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestInput );
    CPPUNIT_TEST( testElement );
    CPPUNIT_TEST_SUITE_END();

  private:
    Input* input_;

  public:
    void setUp() {
        this->input_ = new Input( "tests/labStrains_first100_PLAF.txt",
                                  "tests/PG0390_first100ref.txt",
                                  "tests/PG0390_first100alt.txt",
                                  (size_t)5);
        }

    void tearDown() {
        delete input_;
    }

    void testElement(){
        CPPUNIT_ASSERT_EQUAL((size_t)100, this->input_->nLoci_);
        //"Pf3D7_01_v3"	93157	0.0190612159917058
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0190612159917058, this->input_->plaf[0], 0.000000000001);
        //"Pf3D7_01_v3"	93157	85
        CPPUNIT_ASSERT_EQUAL(85.0, this->input_->refCount[0]);
        //"Pf3D7_01_v3"	93157	0
        CPPUNIT_ASSERT_EQUAL(0.0, this->input_->altCount[0]);

        //"Pf3D7_01_v3"	95518	0.687463394087723
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.687463394087723, this->input_->plaf[4], 0.000000000001);
        //"Pf3D7_01_v3"	95518	156
        CPPUNIT_ASSERT_EQUAL(156.0, this->input_->refCount[4]);
        //"Pf3D7_01_v3"	95518	46
        CPPUNIT_ASSERT_EQUAL(46.0, this->input_->altCount[4]);

        //"Pf3D7_01_v3"	113396	0.0207016179419847
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0207016179419847, this->input_->plaf[20], 0.000000000001);
        //"Pf3D7_01_v3"	113396	177
        CPPUNIT_ASSERT_EQUAL(177.0, this->input_->refCount[20]);
        //"Pf3D7_01_v3"	113396	0
        CPPUNIT_ASSERT_EQUAL(0.0, this->input_->altCount[20]);

        //"Pf3D7_01_v3"	180270	0.426732645350475
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.426732645350475, this->input_->plaf[99], 0.000000000001);
        //"Pf3D7_01_v3"	180270	199
        CPPUNIT_ASSERT_EQUAL(199.0, this->input_->refCount[99]);
        //"Pf3D7_01_v3"	180270	0
        CPPUNIT_ASSERT_EQUAL(0.0, this->input_->altCount[99]);
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestInput );
