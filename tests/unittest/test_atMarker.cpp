#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "atMarker.hpp"

class TestTxtReader : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestTxtReader );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( checkSizeBefore );
    CPPUNIT_TEST( checkInfo );
    CPPUNIT_TEST( checkRemoveMarkers );
    CPPUNIT_TEST( checkSizeAfter );
    CPPUNIT_TEST_SUITE_END();

  private:
    TxtReader * atMarker_;
    InputMarker * altCount_;
    ExcludeMarker* excludedMarkers_;

  public:
    void setUp() {
        this->atMarker_ = new TxtReader();
        this->altCount_ = new InputMarker();
        this->excludedMarkers_ = new ExcludeMarker();
        this->altCount_->readFromFile("tests/testData/atMarkerForTesting.txt" );
        this->excludedMarkers_->readFromFile("tests/testData/excludeFortesting.txt" );
    }

    void tearDown() {
        delete atMarker_;
        delete altCount_;
        delete excludedMarkers_;
    }

    void testMainConstructor(){ }

    void checkSizeBefore(){
        CPPUNIT_ASSERT_EQUAL ( (size_t)100, this->altCount_->info_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)100, this->altCount_->nLoci_ );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->altCount_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->altCount_->position_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)12, this->altCount_->position_[0].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)16, this->altCount_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)17, this->altCount_->position_[2].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)24, this->altCount_->position_[3].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)20, this->altCount_->position_[4].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)11, this->altCount_->position_[5].size() );

        CPPUNIT_ASSERT_EQUAL ( (size_t)0, this->excludedMarkers_->info_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)3, this->excludedMarkers_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)3, this->excludedMarkers_->position_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, this->excludedMarkers_->position_[0].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)1, this->excludedMarkers_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)2, this->excludedMarkers_->position_[2].size() );
    }

    void checkInfo(){
        //Pf3D7_01_v3	93157	0
        CPPUNIT_ASSERT_EQUAL ( (double)93157, this->altCount_->position_[0][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->info_[0] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->content_[0][0] );

        //Pf3D7_01_v3	95518	46
        CPPUNIT_ASSERT_EQUAL ( (double)95518, this->altCount_->position_[0][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)46, this->altCount_->info_[6-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)46, this->altCount_->content_[6-2][0] );

        //Pf3D7_02_v3	100608	35
        CPPUNIT_ASSERT_EQUAL ( (double)100608, this->altCount_->position_[1][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)35, this->altCount_->info_[14-2] );  // it is at line 14 in the file, minus 1 for the header, minus another 1 for 0-indexing
        CPPUNIT_ASSERT_EQUAL ( (double)35, this->altCount_->content_[14-2][0] );

        //Pf3D7_02_v3	107823	0
        CPPUNIT_ASSERT_EQUAL ( (double)107823, this->altCount_->position_[1][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->info_[18-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->content_[18-2][0] );

        //Pf3D7_02_v3	111040	20
        CPPUNIT_ASSERT_EQUAL ( (double)111040, this->altCount_->position_[1][5] );
        CPPUNIT_ASSERT_EQUAL ( (double)20, this->altCount_->info_[19-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)20, this->altCount_->content_[19-2][0] );

        //Pf3D7_02_v3	114473	49
        CPPUNIT_ASSERT_EQUAL ( (double)114473, this->altCount_->position_[1][9] );
        CPPUNIT_ASSERT_EQUAL ( (double)49, this->altCount_->info_[23-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)49, this->altCount_->content_[23-2][0] );

        //Pf3D7_04_v3	144877	47
        CPPUNIT_ASSERT_EQUAL ( (double)144877, this->altCount_->position_[3][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)47, this->altCount_->info_[47-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)47, this->altCount_->content_[47-2][0] );

        //Pf3D7_06_v3	180170	63
        CPPUNIT_ASSERT_EQUAL ( (double)180170, this->altCount_->position_[5][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)63, this->altCount_->info_[91-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)63, this->altCount_->content_[92-2][0] );

        //Pf3D7_06_v3	180183	60
        CPPUNIT_ASSERT_EQUAL ( (double)180183, this->altCount_->position_[5][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)60, this->altCount_->info_[95-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)60, this->altCount_->content_[95-2][0] );
    }

    void checkRemoveMarkers(){
        CPPUNIT_ASSERT_NO_THROW ( this->altCount_->removeMarkers (excludedMarkers_) );
    }

    void checkSizeAfter(){
        this->altCount_->removeMarkers (excludedMarkers_);
        CPPUNIT_ASSERT_EQUAL ( (size_t)93, this->altCount_->info_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)93, this->altCount_->nLoci_ );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->altCount_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->altCount_->position_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)12, this->altCount_->position_[0].size() );
        //Pf3D7_02_v3	100608
        //Pf3D7_02_v3	107823
        //Pf3D7_02_v3	111040
        //Pf3D7_02_v3	114473
        CPPUNIT_ASSERT_EQUAL ( (size_t)(16-4), this->altCount_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)17, this->altCount_->position_[2].size() );
        //Pf3D7_04_v3	144877
        CPPUNIT_ASSERT_EQUAL ( (size_t)(24-1), this->altCount_->position_[3].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)20, this->altCount_->position_[4].size() );
        //Pf3D7_06_v3	180170
        //Pf3D7_06_v3	180183
        CPPUNIT_ASSERT_EQUAL ( (size_t)(11-2), this->altCount_->position_[5].size() );

        CPPUNIT_ASSERT_EQUAL ( (size_t)0, this->excludedMarkers_->info_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)3, this->excludedMarkers_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)3, this->excludedMarkers_->position_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, this->excludedMarkers_->position_[0].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)1, this->excludedMarkers_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)2, this->excludedMarkers_->position_[2].size() );

        //Pf3D7_01_v3	93157	0
        CPPUNIT_ASSERT_EQUAL ( (double)93157, this->altCount_->position_[0][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->info_[0] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->content_[0][0] );

        //Pf3D7_01_v3	95518	46
        CPPUNIT_ASSERT_EQUAL ( (double)95518, this->altCount_->position_[0][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)46, this->altCount_->info_[6-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)46, this->altCount_->content_[6-2][0] );

        //Pf3D7_02_v3	101269	33
        CPPUNIT_ASSERT_EQUAL ( (double)101269, this->altCount_->position_[1][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)33, this->altCount_->info_[14-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)33, this->altCount_->content_[14-2][0] );

        //Pf3D7_02_v3	112038	28
        CPPUNIT_ASSERT_EQUAL ( (double)112038, this->altCount_->position_[1][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)28, this->altCount_->info_[18-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)28, this->altCount_->content_[18-2][0] );

        //Pf3D7_02_v3	113396	0
        CPPUNIT_ASSERT_EQUAL ( (double)113396, this->altCount_->position_[1][5] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->info_[19-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->content_[19-2][0] );

        //Pf3D7_02_v3	121693	0
        CPPUNIT_ASSERT_EQUAL ( (double)121693, this->altCount_->position_[1][9] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->info_[23-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->content_[23-2][0] );

        //Pf3D7_04_v3	149598	0
        CPPUNIT_ASSERT_EQUAL ( (double)149598, this->altCount_->position_[3][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->info_[47-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->content_[47-2][0] );

        //Pf3D7_06_v3	180217	0
        CPPUNIT_ASSERT_EQUAL ( (double)180217, this->altCount_->position_[5][5] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->info_[91-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->content_[91-2][0] );

        //Pf3D7_06_v3	180270	0
        CPPUNIT_ASSERT_EQUAL ( (double)180270, this->altCount_->position_[5][8] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->info_[94-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->altCount_->content_[94-2][0] );
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestTxtReader );
