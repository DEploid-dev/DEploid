#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../../src/panel.hpp"

class TestPanel : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestPanel );
    CPPUNIT_TEST( testElement );
    CPPUNIT_TEST_SUITE_END();

  private:
    Panel* panel_;

  public:
    void setUp() {
        this->panel_ = new Panel("tests/lab_first100_Panel.txt");
    }

    void tearDown() {
        delete panel_;
    }

    void testElement(){
        //Pf3D7_01_v3,93157,0,0,0,1
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[0][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[0][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[0][2], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[0][3], 1.0);

        //Pf3D7_01_v3,95518,0,1,1,0
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[4][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[4][1], 1.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[4][2], 1.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[4][3], 0.0);

        //Pf3D7_01_v3,113396,0,0,0,1
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[20][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[20][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[20][2], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[20][3], 1.0);

        //Pf3D7_01_v3,180270,0,0,0,1
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[99][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[99][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[99][2], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content[99][3], 1.0);
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestPanel );
