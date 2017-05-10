#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "panel.hpp"

class TestPanel : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestPanel );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( checkContent );
    CPPUNIT_TEST( checkChrom );
    CPPUNIT_TEST( checkChromIndex );
    CPPUNIT_TEST( checkPOS );
    CPPUNIT_TEST( checkRecombProb );
    CPPUNIT_TEST( checkRecombProbEach );
    CPPUNIT_TEST( checkRecRec );
    CPPUNIT_TEST( checkRecNoRec );
    CPPUNIT_TEST( checkNoRecNoRec );
    CPPUNIT_TEST( checkForbid );
    CPPUNIT_TEST( checkRecomb1 );
    CPPUNIT_TEST( checkRecomb0 );
    CPPUNIT_TEST( checkForbidRecomb1 );
    CPPUNIT_TEST( checkForbidRecomb0 );
    CPPUNIT_TEST( testLociNumberUnequal );
    CPPUNIT_TEST( checkUpdatingReferencePanel);
    CPPUNIT_TEST( checkExamplePanel );
    CPPUNIT_TEST_SUITE_END();

  private:
    Panel* panel1_;
    Panel* panel2_;
    Panel* panel3_;
    Panel* panel4_;

    string panelName_;
    double epsilon3;
  public:
    void setUp() {
        this->epsilon3 = 0.000000000001;
        this->panelName_ = "data/testData/testingPanel.txt";
        // in R: panel = read.csv("clonalPanel.csv", header = T)
        this->panel1_ = new Panel();
        this->panel1_->readFromFile(panelName_.c_str() );
        this->panel1_->computeRecombProbs( 15000.0, 20.0, false, 0, false); // forbid copy from same = false by default!

        this->panel2_ = new Panel();
        this->panel2_->readFromFile( panelName_.c_str() );

        this->panel3_ = new Panel();
        this->panel3_->buildExamplePanel1();

        this->panel4_ = new Panel();
        this->panel4_->buildExamplePanel2();
    }


    void tearDown() {
        delete panel1_;
        delete panel2_;
        delete panel3_;
        delete panel4_;
    }


    void testMainConstructor(){ }


    void checkContent(){
        //Pf3D7_01_v3,93157,0,0,0,1
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[0][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[0][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[0][2], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[0][3], 1.0);

        //Pf3D7_01_v3,95518,0,1,1,0
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[4][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[4][1], 1.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[4][2], 1.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[4][3], 0.0);

        //Pf3D7_01_v3,113396,0,0,0,1
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[20][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[20][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[20][2], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[20][3], 1.0);

        //Pf3D7_01_v3,180270,0,0,0,1
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[99][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[99][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[99][2], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[99][3], 1.0);

        //Pf3D7_09_v3,89744,0,0,0,1
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[8098][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[8098][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[8098][2], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->panel1_->content_[8098][3], 1.0);
    }


    void checkChromIndex(){
        CPPUNIT_ASSERT_EQUAL( (int)13, this->panel1_->tmpChromInex_ );
    }


    void checkChrom(){
        CPPUNIT_ASSERT_EQUAL( (size_t)14, this->panel1_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_01_v3"), this->panel1_->chrom_[0] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_02_v3"), this->panel1_->chrom_[1] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_03_v3"), this->panel1_->chrom_[2] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_10_v3"), this->panel1_->chrom_[9] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_11_v3"), this->panel1_->chrom_[10] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_14_v3"), this->panel1_->chrom_.back() );

        CPPUNIT_ASSERT_EQUAL( (size_t)14, this->panel2_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_01_v3"), this->panel2_->chrom_[0] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_02_v3"), this->panel2_->chrom_[1] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_03_v3"), this->panel2_->chrom_[2] );

        CPPUNIT_ASSERT_EQUAL( (size_t)1, this->panel3_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_01_v3"), this->panel3_->chrom_[0] );

        CPPUNIT_ASSERT_EQUAL( (size_t)3, this->panel4_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_01_v3"), this->panel4_->chrom_[0] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_02_v3"), this->panel4_->chrom_[1] );
        CPPUNIT_ASSERT_EQUAL( string("Pf3D7_03_v3"), this->panel4_->chrom_[2] );
    }


    void checkPOS(){
        CPPUNIT_ASSERT_EQUAL( (size_t)14, this->panel1_->position_.size() );

        CPPUNIT_ASSERT_EQUAL( (size_t)590, this->panel1_->position_[0].size() );
        CPPUNIT_ASSERT_EQUAL( (int)93157, this->panel1_->position_[0][0] );
        CPPUNIT_ASSERT_EQUAL( (int)95518, this->panel1_->position_[0][4] );
        CPPUNIT_ASSERT_EQUAL( (int)113396, this->panel1_->position_[0][20] );
        CPPUNIT_ASSERT_EQUAL( (int)180270, this->panel1_->position_[0][99] );

        CPPUNIT_ASSERT_EQUAL( (size_t)705, this->panel1_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)742, this->panel1_->position_[2].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1547, this->panel1_->position_[3].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)937, this->panel1_->position_[4].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)956, this->panel1_->position_[5].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1432, this->panel1_->position_[6].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1180, this->panel1_->position_[7].size() );

        CPPUNIT_ASSERT_EQUAL( (size_t)1079, this->panel1_->position_[8].size() );
        // in R: panel[panel$CHROM=="Pf3D7_09_v3",][10,]
        CPPUNIT_ASSERT_EQUAL( (int)89744, this->panel1_->position_[8][9] );

        CPPUNIT_ASSERT_EQUAL( (size_t)1275, this->panel1_->position_[9].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1385, this->panel1_->position_[10].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1240, this->panel1_->position_[11].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)1844, this->panel1_->position_[12].size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)2203, this->panel1_->position_[13].size() );

        CPPUNIT_ASSERT_EQUAL( (size_t)14, this->panel2_->position_.size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)590, this->panel2_->position_[0].size() );
        CPPUNIT_ASSERT_EQUAL( (int)93157, this->panel2_->position_[0][0] );
        CPPUNIT_ASSERT_EQUAL( (int)95518, this->panel2_->position_[0][4] );
        CPPUNIT_ASSERT_EQUAL( (int)113396, this->panel2_->position_[0][20] );
        CPPUNIT_ASSERT_EQUAL( (int)180270, this->panel2_->position_[0][99] );

        CPPUNIT_ASSERT_EQUAL( (size_t)705, this->panel2_->position_[1].size() );

        CPPUNIT_ASSERT_EQUAL( (size_t)1, this->panel3_->position_.size() );
        CPPUNIT_ASSERT_EQUAL( (size_t)3, this->panel4_->position_.size() );

    }


    void checkRecombProb(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029, this->panel1_->pRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042, this->panel1_->pRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074, this->panel1_->pRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115, this->panel1_->pRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896, this->panel1_->pRec_[4], epsilon3);

        // in R: cumsum(table(panel$CHROM))-1
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[589]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[1294]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[2036]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[3583]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[4520]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[5476]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[6908]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[8088]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[9167]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[10442]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[11827]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[13067]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[14911]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel1_->pRec_[17114]);
    }


    void checkRecombProbEach(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029/4.0, this->panel1_->pRecEachHap_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0, this->panel1_->pRecEachHap_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0, this->panel1_->pRecEachHap_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0, this->panel1_->pRecEachHap_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0, this->panel1_->pRecEachHap_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[589], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[1294], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[2036], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[3583], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[4520], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[5476], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[6908], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[8088], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[9167], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[10442], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[11827], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[13067], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[14911], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel1_->pRecEachHap_[17114], epsilon3);
    }


    void checkRecRec(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029/4.0*0.016725220801029/4.0, this->panel1_->pRecRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0*0.000493211664453042/4.0, this->panel1_->pRecRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0*0.000373263653116074/4.0, this->panel1_->pRecRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0*0.0136526127138115/4.0, this->panel1_->pRecRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0*0.00151884538507896/4.0, this->panel1_->pRecRec_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[589], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[1294], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[2036], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[3583], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[4520], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[5476], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[6908], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[8088], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[9167], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[10442], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[11827], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[13067], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[14911], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel1_->pRecRec_[17114], epsilon3);
    }


    void checkRecNoRec(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029/4.0*(1-0.016725220801029), this->panel1_->pRecNoRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0*(1-0.000493211664453042), this->panel1_->pRecNoRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0*(1-0.000373263653116074), this->panel1_->pRecNoRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0*(1-0.0136526127138115), this->panel1_->pRecNoRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0*(1-0.00151884538507896), this->panel1_->pRecNoRec_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[589], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[1294], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[2036], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[3583], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[4520], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[5476], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[6908], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[8088], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[9167], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[10442], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[11827], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[13067], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[14911], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pRecNoRec_[17114], epsilon3);
    }


    void checkNoRecNoRec(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.016725220801029)*(1-0.016725220801029), this->panel1_->pNoRecNoRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.000493211664453042)*(1-0.000493211664453042), this->panel1_->pNoRecNoRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.000373263653116074)*(1-0.000373263653116074), this->panel1_->pNoRecNoRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.0136526127138115)*(1-0.0136526127138115), this->panel1_->pNoRecNoRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.00151884538507896)*(1-0.00151884538507896), this->panel1_->pNoRecNoRec_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[589], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[1294], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[2036], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[3583], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[4520], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[5476], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[6908], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[8088], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[9167], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[10442], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[11827], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[13067], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[14911], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel1_->pNoRecNoRec_[17114], epsilon3);
    }


    void testLociNumberUnequal(){
        CPPUNIT_ASSERT_THROW ( this->panel1_->checkForExceptions(17114, this->panelName_), LociNumberUnequal );
        CPPUNIT_ASSERT_THROW ( this->panel1_->checkForExceptions(100, this->panelName_), LociNumberUnequal );
        CPPUNIT_ASSERT_NO_THROW ( this->panel1_->checkForExceptions(17115, this->panelName_) );

        // Test if recombination probabilitiy were not computed first, throw exception
        CPPUNIT_ASSERT_THROW ( this->panel2_->checkForExceptions(17115, this->panelName_), LociNumberUnequal );
        this->panel2_->computeRecombProbs( 15000.0, 10.0, false, 0, true); // forbid copy from same = false by default!

        // After recombination probability was computed, no error thrown
        CPPUNIT_ASSERT_NO_THROW ( this->panel2_->checkForExceptions(17115, this->panelName_) );
    }


    void checkForbid(){
        this->panel2_->computeRecombProbs( 15000.0, 20.0, false, 0, true); // forbid copy from same = false by default!

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029/4.0, this->panel2_->pRecEachHap_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0, this->panel2_->pRecEachHap_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0, this->panel2_->pRecEachHap_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0, this->panel2_->pRecEachHap_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0, this->panel2_->pRecEachHap_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029/4.0*0.016725220801029/3.0, this->panel2_->pRecRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0*0.000493211664453042/3.0, this->panel2_->pRecRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0*0.000373263653116074/3.0, this->panel2_->pRecRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0*0.0136526127138115/3.0, this->panel2_->pRecRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0*0.00151884538507896/3.0, this->panel2_->pRecRec_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[589], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[1294], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[2036], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[3583], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[4520], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[5476], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[6908], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[8088], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[9167], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[10442], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[11827], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[13067], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[14911], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3, this->panel2_->pRecRec_[17114], epsilon3);
    }


    void checkRecomb1(){
        this->panel2_->computeRecombProbs( 15000.0, 20.0, true, 1, false); // forbid copy from same = false by default!

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel2_->pRecEachHap_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel2_->pRecEachHap_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel2_->pRecEachHap_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel2_->pRecEachHap_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel2_->pRecEachHap_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[589], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[1294], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[2036], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[3583], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[4520], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[5476], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[6908], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[8088], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[9167], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[10442], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[11827], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[13067], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[14911], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[17114], epsilon3);
    }


    void checkRecomb0(){
        this->panel2_->computeRecombProbs( 15000.0, 10.0, true, 0, false); // forbid copy from same = false by default!
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecEachHap_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecEachHap_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecEachHap_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecEachHap_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecEachHap_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecRec_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[589], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[1294], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[2036], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[3583], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[4520], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[5476], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[6908], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[8088], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[9167], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[10442], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[11827], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[13067], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[14911], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/4.0, this->panel2_->pRecRec_[17114], epsilon3);
    }


    void checkForbidRecomb1(){
        this->panel2_->computeRecombProbs( 15000.0, 10.0, true, 1, true); // forbid copy from same = false by default!

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel2_->pRecEachHap_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel2_->pRecEachHap_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel2_->pRecEachHap_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel2_->pRecEachHap_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel2_->pRecEachHap_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[589], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[1294], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[2036], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[3583], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[4520], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[5476], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[6908], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[8088], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[9167], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[10442], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[11827], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[13067], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[14911], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[17114], epsilon3);
    }


    void checkForbidRecomb0(){
        this->panel2_->computeRecombProbs( 15000.0, 10.0, true, 0.0, true); // forbid copy from same = false by default!

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecEachHap_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecEachHap_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecEachHap_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecEachHap_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecEachHap_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, this->panel2_->pRecRec_[4], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[589], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[1294], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[2036], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[3583], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[4520], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[5476], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[6908], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[8088], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[9167], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[10442], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[11827], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[13067], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[14911], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25/3.0, this->panel2_->pRecRec_[17114], epsilon3);
    }

    void checkExamplePanel(){
        CPPUNIT_ASSERT_NO_THROW ( this->panel3_->computeRecombProbs( 15000.0, 20.0, false, 0, false) ); // forbid copy from same = false by default!
        CPPUNIT_ASSERT_NO_THROW ( this->panel3_->checkForExceptions(7, "Panel 3") );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029, this->panel3_->pRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042, this->panel3_->pRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074, this->panel3_->pRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115, this->panel3_->pRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896, this->panel3_->pRec_[4], epsilon3);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel3_->pRec_[6]);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029/4.0, this->panel3_->pRecEachHap_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0, this->panel3_->pRecEachHap_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0, this->panel3_->pRecEachHap_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0, this->panel3_->pRecEachHap_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0, this->panel3_->pRecEachHap_[4], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel3_->pRecEachHap_[6], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029/4.0*0.016725220801029/4.0, this->panel3_->pRecRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0*0.000493211664453042/4.0, this->panel3_->pRecRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0*0.000373263653116074/4.0, this->panel3_->pRecRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0*0.0136526127138115/4.0, this->panel3_->pRecRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0*0.00151884538507896/4.0, this->panel3_->pRecRec_[4], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel3_->pRecRec_[6], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.016725220801029/4.0*(1-0.016725220801029), this->panel3_->pRecNoRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0*(1-0.000493211664453042), this->panel3_->pRecNoRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0*(1-0.000373263653116074), this->panel3_->pRecNoRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0*(1-0.0136526127138115), this->panel3_->pRecNoRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0*(1-0.00151884538507896), this->panel3_->pRecNoRec_[4], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel3_->pRecNoRec_[6], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.016725220801029)*(1-0.016725220801029), this->panel3_->pNoRecNoRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.000493211664453042)*(1-0.000493211664453042), this->panel3_->pNoRecNoRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.000373263653116074)*(1-0.000373263653116074), this->panel3_->pNoRecNoRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.0136526127138115)*(1-0.0136526127138115), this->panel3_->pNoRecNoRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.00151884538507896)*(1-0.00151884538507896), this->panel3_->pNoRecNoRec_[4], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel3_->pNoRecNoRec_[6], epsilon3);

        CPPUNIT_ASSERT_NO_THROW ( this->panel4_->computeRecombProbs( 15000.0, 20.0, false, 0, false) ); // forbid copy from same = false by default!
        CPPUNIT_ASSERT_NO_THROW ( this->panel4_->checkForExceptions(7, "Panel 4") );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, this->panel4_->pRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042, this->panel4_->pRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074, this->panel4_->pRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115, this->panel4_->pRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896, this->panel4_->pRec_[4], epsilon3);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel4_->pRec_[5]);
        CPPUNIT_ASSERT_EQUAL(1.0, this->panel4_->pRec_[6]);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel4_->pRecEachHap_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0, this->panel4_->pRecEachHap_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0, this->panel4_->pRecEachHap_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0, this->panel4_->pRecEachHap_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0, this->panel4_->pRecEachHap_[4], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel4_->pRecEachHap_[5], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, this->panel4_->pRecEachHap_[6], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel4_->pRecRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0*0.000493211664453042/4.0, this->panel4_->pRecRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0*0.000373263653116074/4.0, this->panel4_->pRecRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0*0.0136526127138115/4.0, this->panel4_->pRecRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0*0.00151884538507896/4.0, this->panel4_->pRecRec_[4], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel4_->pRecRec_[5], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*0.25, this->panel4_->pRecRec_[6], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel4_->pRecNoRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000493211664453042/4.0*(1-0.000493211664453042), this->panel4_->pRecNoRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.000373263653116074/4.0*(1-0.000373263653116074), this->panel4_->pRecNoRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0136526127138115/4.0*(1-0.0136526127138115), this->panel4_->pRecNoRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00151884538507896/4.0*(1-0.00151884538507896), this->panel4_->pRecNoRec_[4], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel4_->pRecNoRec_[5], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel4_->pRecNoRec_[6], epsilon3);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel4_->pNoRecNoRec_[0], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.000493211664453042)*(1-0.000493211664453042), this->panel4_->pNoRecNoRec_[1], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.000373263653116074)*(1-0.000373263653116074), this->panel4_->pNoRecNoRec_[2], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.0136526127138115)*(1-0.0136526127138115), this->panel4_->pNoRecNoRec_[3], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL((1-0.00151884538507896)*(1-0.00151884538507896), this->panel4_->pNoRecNoRec_[4], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel4_->pNoRecNoRec_[5], epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, this->panel4_->pNoRecNoRec_[6], epsilon3);


    }

    void checkUpdatingReferencePanel(){
        vector < vector<double> >  hapForTesting;
        hapForTesting.push_back(vector <double> ({0,1,2,4,8}));
        hapForTesting.push_back(vector <double> ({0,1,2,4,8}));
        hapForTesting.push_back(vector <double> ({0,1,2,4,8}));
        hapForTesting.push_back(vector <double> ({0,1,2,4,8}));
        hapForTesting.push_back(vector <double> ({0,1,2,4,8}));
        hapForTesting.push_back(vector <double> ({0,1,2,4,8}));
        hapForTesting.push_back(vector <double> ({0,1,2,4,8}));

        CPPUNIT_ASSERT_NO_THROW(this->panel3_->initializeUpdatePanel(8));
        for ( size_t i = 0; i < 7; i++){
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i].size(), (size_t)8);
        }

        CPPUNIT_ASSERT_NO_THROW(this->panel3_->updatePanelWithHaps(8, 1, hapForTesting));
        for ( size_t i = 0; i < 7; i++){
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][4], (double)0);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][5], (double)2);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][6], (double)4);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][7], (double)8);
        }

        CPPUNIT_ASSERT_NO_THROW(this->panel3_->updatePanelWithHaps(8, 4, hapForTesting));
        for ( size_t i = 0; i < 7; i++){
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][4], (double)0);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][5], (double)1);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][6], (double)2);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][7], (double)4);
        }

        CPPUNIT_ASSERT_NO_THROW(this->panel3_->updatePanelWithHaps(8, 2, hapForTesting));
        for ( size_t i = 0; i < 7; i++){
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][4], (double)0);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][5], (double)1);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][6], (double)4);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][7], (double)8);
        }

        CPPUNIT_ASSERT_NO_THROW(this->panel3_->updatePanelWithHaps(8, 0, hapForTesting));
        for ( size_t i = 0; i < 7; i++){
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][4], (double)1);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][5], (double)2);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][6], (double)4);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][7], (double)8);
        }

        CPPUNIT_ASSERT_NO_THROW(this->panel3_->updatePanelWithHaps(8, 3, hapForTesting));
        for ( size_t i = 0; i < 7; i++){
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][4], (double)0);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][5], (double)1);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][6], (double)2);
            CPPUNIT_ASSERT_EQUAL(this->panel3_->content_[i][7], (double)8);
        }

    }

};


class TestInitialHaplotypes : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestInitialHaplotypes );
    CPPUNIT_TEST( testMainConstructor );
    //CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST_SUITE_END();

  private:
    InitialHaplotypes* hap_;

    string hapName_;
  public:
    void setUp() {
        this->hapName_ = "data/testData/PG0390-C.test.nopanel.hap";
        this->hap_ = new InitialHaplotypes();
        this->hap_->readFromFile(hapName_.c_str() );
    }

    void tearDown() {
        delete hap_;
    }
  public:
    void testMainConstructor(){
        CPPUNIT_ASSERT_EQUAL(this->hap_->chrom_.size(), (size_t)14);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_.size(), (size_t)14);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_.size(), (size_t)14);

        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[0], (size_t)0);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[1], (size_t)204);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[2], (size_t)216);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[3], (size_t)235);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[4], (size_t)275);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[5], (size_t)299);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[6], (size_t)320);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[7], (size_t)362);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[8], (size_t)393);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[9], (size_t)429);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[10], (size_t)450);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[11], (size_t)478);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[12], (size_t)504);
        CPPUNIT_ASSERT_EQUAL(this->hap_->indexOfChromStarts_[13], (size_t)545);

        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[0].size(), (size_t)204);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[0][0], (int)93157);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[0][1], (int)94422);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[0][2], (int)94459);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[1].size(), (size_t)12);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[2].size(), (size_t)19);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[3].size(), (size_t)40);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[4].size(), (size_t)24);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[5].size(), (size_t)21);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[5][0], (int)120625);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[5][1], (int)183045);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[5][2], (int)244644);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[6].size(), (size_t)42);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[7].size(), (size_t)31);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[8].size(), (size_t)36);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[9].size(), (size_t)21);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[10].size(), (size_t)28);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[11].size(), (size_t)26);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[12].size(), (size_t)41);
        CPPUNIT_ASSERT_EQUAL(this->hap_->position_[13].size(), (size_t)49);

        //head -5 data/testData/PG0390-C.test.nopanel.hap | tail -n1
        //Pf3D7_01_v3	94487	0	0	1	0	0

        CPPUNIT_ASSERT_EQUAL(this->hap_->content_.size(), (size_t)594);
        CPPUNIT_ASSERT_EQUAL(this->hap_->content_[3][0], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->hap_->content_[3][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->hap_->content_[3][2], 1.0);
        CPPUNIT_ASSERT_EQUAL(this->hap_->content_[3][3], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->hap_->content_[3][4], 0.0);

        //head -195 data/testData/PG0390-C.test.nopanel.hap | tail -n1
        //Pf3D7_01_v3	319673	1	0	1	0	0
        CPPUNIT_ASSERT_EQUAL(this->hap_->content_[193][0], 1.0);
        CPPUNIT_ASSERT_EQUAL(this->hap_->content_[193][1], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->hap_->content_[193][2], 1.0);
        CPPUNIT_ASSERT_EQUAL(this->hap_->content_[193][3], 0.0);
        CPPUNIT_ASSERT_EQUAL(this->hap_->content_[193][4], 0.0);
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestPanel );
CPPUNIT_TEST_SUITE_REGISTRATION( TestInitialHaplotypes );
