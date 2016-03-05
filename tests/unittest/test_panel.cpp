#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../../src/panel.hpp"

class TestPanel : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestPanel );
    CPPUNIT_TEST( checkContent );
    CPPUNIT_TEST( checkChrom );
    CPPUNIT_TEST( checkChromIndex );
    CPPUNIT_TEST( checkPOS );
    CPPUNIT_TEST( checkPOS2 );
    CPPUNIT_TEST( checkRecombProb );
    CPPUNIT_TEST( checkRecombProbEach );
    CPPUNIT_TEST_SUITE_END();

  private:
    Panel* panel_;

  public:
    void setUp() {
        // in R: panel = read.csv("clonalPanel.csv", header = T)
        this->panel_ = new Panel("labStrains/clonalPanel.csv");
    }

    void tearDown() {
        delete panel_;
    }

    void checkContent(){
        //Pf3D7_01_v3,93157,0,0,0,1
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[0][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[0][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[0][2], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[0][3], 1.0);

        //Pf3D7_01_v3,95518,0,1,1,0
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[4][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[4][1], 1.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[4][2], 1.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[4][3], 0.0);

        //Pf3D7_01_v3,113396,0,0,0,1
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[20][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[20][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[20][2], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[20][3], 1.0);

        //Pf3D7_01_v3,180270,0,0,0,1
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[99][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[99][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[99][2], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[99][3], 1.0);

        //Pf3D7_09_v3,89744,0,0,0,1
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[8098][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[8098][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[8098][2], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel_->content_[8098][3], 1.0);
    }

    void checkChromIndex(){
        CPPUNIT_ASSERT_EQUAL( (int)13, this->panel_->tmpChromInex_ );
    }

    void checkChrom(){
        CPPUNIT_ASSERT_EQUAL( (size_t)14, this->panel_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_01_v3"), this->panel_->chrom_[0] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_02_v3"), this->panel_->chrom_[1] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_03_v3"), this->panel_->chrom_[2] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_10_v3"), this->panel_->chrom_[9] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_11_v3"), this->panel_->chrom_[10] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_14_v3"), this->panel_->chrom_.back() );
    }

    void checkPOS(){
        CPPUNIT_ASSERT_EQUAL( (size_t)14, this->panel_->position_.size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)590, this->panel_->position_[0].size() );
        CPPUNIT_ASSERT_EQUAL( (double)93157, this->panel_->position_[0][0] );
        CPPUNIT_ASSERT_EQUAL( (double)95518, this->panel_->position_[0][4] );
        CPPUNIT_ASSERT_EQUAL( (double)113396, this->panel_->position_[0][20] );
        CPPUNIT_ASSERT_EQUAL( (double)180270, this->panel_->position_[0][99] );
        CPPUNIT_ASSERT_EQUAL( (size_t)705, this->panel_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)742, this->panel_->position_[2].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1547, this->panel_->position_[3].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)937, this->panel_->position_[4].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)956, this->panel_->position_[5].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1432, this->panel_->position_[6].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1180, this->panel_->position_[7].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1079, this->panel_->position_[8].size() );
        // in R: panel[panel$CHROM=="Pf3D7_09_v3",][10,]
        CPPUNIT_ASSERT_EQUAL( (double)89744, this->panel_->position_[8][9] );
        CPPUNIT_ASSERT_EQUAL( (size_t)1275, this->panel_->position_[9].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1385, this->panel_->position_[10].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1240, this->panel_->position_[11].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1844, this->panel_->position_[12].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)2203, this->panel_->position_[13].size() );
    }

    void checkPOS2(){
        CPPUNIT_ASSERT_EQUAL( (size_t)14, this->panel_->position_.size() );
    }

    void checkRecombProb(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029, this->panel_->pRec_[0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042, this->panel_->pRec_[1], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074, this->panel_->pRec_[2], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115, this->panel_->pRec_[3], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896, this->panel_->pRec_[4], 0.000000000001);

        // in R: cumsum(table(panel$CHROM))-1
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[589]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[1294]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[2036]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[3583]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[4520]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[5476]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[6908]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[8088]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[9167]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[10442]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[11827]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[13067]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[14911]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel_->pRec_[17114]);
    }

    void checkRecombProbEach(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029/4.0, this->panel_->pRecEachHap_[0], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0, this->panel_->pRecEachHap_[1], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0, this->panel_->pRecEachHap_[2], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0, this->panel_->pRecEachHap_[3], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0, this->panel_->pRecEachHap_[4], 0.000000000001);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[589], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[1294], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[2036], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[3583], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[4520], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[5476], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[6908], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[8088], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[9167], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[10442], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[11827], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[13067], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[14911], 0.000000000001);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel_->pRecEachHap_[17114], 0.000000000001);
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestPanel );
