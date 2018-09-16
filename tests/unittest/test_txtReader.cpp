#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "txtReader.hpp"

class TestTxtReader : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestTxtReader );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( checkSizeBefore );
    CPPUNIT_TEST( checkInfo );
    CPPUNIT_TEST( checkRemoveMarkers );
    CPPUNIT_TEST( checkSizeAfter );
    CPPUNIT_TEST( checkSortedPositions );
    CPPUNIT_TEST( checkBadConversion );
    CPPUNIT_TEST( checkBadScientificNotation );
    CPPUNIT_TEST_SUITE_END();

  private:
    TxtReader * txtReader_;
    TxtReader * afterExclude_;
    ExcludeMarker* excludedMarkers_;

  public:
    void setUp() {
        this->txtReader_ = new TxtReader();
        this->afterExclude_ = new TxtReader();
        this->excludedMarkers_ = new ExcludeMarker();
        this->txtReader_->readFromFile("data/testData/txtReaderForTesting.txt" );
        this->afterExclude_->readFromFile("data/testData/txtReaderForTestingAfterExclude.txt");
        this->excludedMarkers_->readFromFile("data/testData/txtReaderForTestingToBeExclude.txt" );
    }

    void tearDown() {
        delete txtReader_;
        delete afterExclude_;
        delete excludedMarkers_;
    }

    void testMainConstructor(){ }

    void checkBadConversion(){
        TxtReader tmp;
        CPPUNIT_ASSERT_THROW ( tmp.readFromFile("data/testData/bad.plaf.txt"), BadConversion );
    }

    void checkBadScientificNotation(){
        TxtReader tmp2;
        CPPUNIT_ASSERT_THROW ( tmp2.readFromFile("data/testData/bad.plaf_scientific.txt"), BadScientificNotation );

        TxtReader tmp3;
        CPPUNIT_ASSERT_THROW ( tmp3.readFromFile("data/testData/bad.plaf_scientificE.txt"), BadScientificNotation );
    }

    void checkSortedPositions() {
        TxtReader tmp;
        CPPUNIT_ASSERT_THROW ( tmp.readFromFile("data/testData/bad.plaf_badpos.txt"), PositionUnsorted );
    }

    void checkSizeBefore(){
        CPPUNIT_ASSERT_EQUAL ( (size_t)100, this->txtReader_->info_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)100, this->txtReader_->nLoci_ );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->txtReader_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->txtReader_->position_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)12, this->txtReader_->position_[0].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)16, this->txtReader_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)17, this->txtReader_->position_[2].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)24, this->txtReader_->position_[3].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)20, this->txtReader_->position_[4].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)11, this->txtReader_->position_[5].size() );

        CPPUNIT_ASSERT_EQUAL ( (size_t)0, this->excludedMarkers_->info_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)3, this->excludedMarkers_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)3, this->excludedMarkers_->position_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, this->excludedMarkers_->position_[0].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)1, this->excludedMarkers_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)2, this->excludedMarkers_->position_[2].size() );

        CPPUNIT_ASSERT_EQUAL ( (size_t)93, this->afterExclude_->info_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)93, this->afterExclude_->nLoci_ );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->afterExclude_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->afterExclude_->position_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)12, this->afterExclude_->position_[0].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)12, this->afterExclude_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)17, this->afterExclude_->position_[2].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)23, this->afterExclude_->position_[3].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)20, this->afterExclude_->position_[4].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)9, this->afterExclude_->position_[5].size() );
    }

    void checkInfo(){
        //Pf3D7_01_v3	93157	0
        CPPUNIT_ASSERT_EQUAL ( (int)93157, this->txtReader_->position_[0][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->info_[0] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->content_[0][0] );

        //Pf3D7_01_v3	95518	46
        CPPUNIT_ASSERT_EQUAL ( (int)95518, this->txtReader_->position_[0][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)46, this->txtReader_->info_[6-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)46, this->txtReader_->content_[6-2][0] );

        //Pf3D7_02_v3	100608	35
        CPPUNIT_ASSERT_EQUAL ( (int)100608, this->txtReader_->position_[1][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)35, this->txtReader_->info_[14-2] );  // it is at line 14 in the file, minus 1 for the header, minus another 1 for 0-indexing
        CPPUNIT_ASSERT_EQUAL ( (double)35, this->txtReader_->content_[14-2][0] );

        //Pf3D7_02_v3	107823	0
        CPPUNIT_ASSERT_EQUAL ( (int)107823, this->txtReader_->position_[1][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->info_[18-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->content_[18-2][0] );

        //Pf3D7_02_v3	111040	20
        CPPUNIT_ASSERT_EQUAL ( (int)111040, this->txtReader_->position_[1][5] );
        CPPUNIT_ASSERT_EQUAL ( (double)20, this->txtReader_->info_[19-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)20, this->txtReader_->content_[19-2][0] );

        //Pf3D7_02_v3	114473	49
        CPPUNIT_ASSERT_EQUAL ( (int)114473, this->txtReader_->position_[1][9] );
        CPPUNIT_ASSERT_EQUAL ( (double)49, this->txtReader_->info_[23-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)49, this->txtReader_->content_[23-2][0] );

        //Pf3D7_04_v3	144877	47
        CPPUNIT_ASSERT_EQUAL ( (int)144877, this->txtReader_->position_[3][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)47, this->txtReader_->info_[47-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)47, this->txtReader_->content_[47-2][0] );

        //Pf3D7_06_v3	180170	63
        CPPUNIT_ASSERT_EQUAL ( (int)180170, this->txtReader_->position_[5][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)63, this->txtReader_->info_[91-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)63, this->txtReader_->content_[92-2][0] );

        //Pf3D7_06_v3	180183	60
        CPPUNIT_ASSERT_EQUAL ( (int)180183, this->txtReader_->position_[5][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)60, this->txtReader_->info_[95-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)60, this->txtReader_->content_[95-2][0] );
    }

    void checkRemoveMarkers(){
        CPPUNIT_ASSERT_NO_THROW ( this->txtReader_->findAndKeepMarkers (excludedMarkers_) );
        CPPUNIT_ASSERT_NO_THROW ( this->afterExclude_->findAndKeepMarkers (excludedMarkers_) );
    }

    void checkSizeAfter(){
        this->txtReader_->findAndKeepMarkers (excludedMarkers_);
        CPPUNIT_ASSERT_EQUAL ( (size_t)93, this->txtReader_->info_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)93, this->txtReader_->nLoci_ );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->txtReader_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->txtReader_->position_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)12, this->txtReader_->position_[0].size() );
        //Pf3D7_02_v3	100608
        //Pf3D7_02_v3	107823
        //Pf3D7_02_v3	111040
        //Pf3D7_02_v3	114473
        CPPUNIT_ASSERT_EQUAL ( (size_t)(16-4), this->txtReader_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)17, this->txtReader_->position_[2].size() );
        //Pf3D7_04_v3	144877
        CPPUNIT_ASSERT_EQUAL ( (size_t)(24-1), this->txtReader_->position_[3].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)20, this->txtReader_->position_[4].size() );
        //Pf3D7_06_v3	180170
        //Pf3D7_06_v3	180183
        CPPUNIT_ASSERT_EQUAL ( (size_t)(11-2), this->txtReader_->position_[5].size() );

        CPPUNIT_ASSERT_EQUAL ( (size_t)0, this->excludedMarkers_->info_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)3, this->excludedMarkers_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)3, this->excludedMarkers_->position_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)4, this->excludedMarkers_->position_[0].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)1, this->excludedMarkers_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)2, this->excludedMarkers_->position_[2].size() );

        //Pf3D7_01_v3	93157	0
        CPPUNIT_ASSERT_EQUAL ( (int)93157, this->txtReader_->position_[0][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->info_[0] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->content_[0][0] );

        //Pf3D7_01_v3	95518	46
        CPPUNIT_ASSERT_EQUAL ( (int)95518, this->txtReader_->position_[0][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)46, this->txtReader_->info_[6-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)46, this->txtReader_->content_[6-2][0] );

        //Pf3D7_02_v3	101269	33
        CPPUNIT_ASSERT_EQUAL ( (int)101269, this->txtReader_->position_[1][0] );
        CPPUNIT_ASSERT_EQUAL ( (double)33, this->txtReader_->info_[14-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)33, this->txtReader_->content_[14-2][0] );

        //Pf3D7_02_v3	112038	28
        CPPUNIT_ASSERT_EQUAL ( (int)112038, this->txtReader_->position_[1][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)28, this->txtReader_->info_[18-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)28, this->txtReader_->content_[18-2][0] );

        //Pf3D7_02_v3	113396	0
        CPPUNIT_ASSERT_EQUAL ( (int)113396, this->txtReader_->position_[1][5] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->info_[19-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->content_[19-2][0] );

        //Pf3D7_02_v3	121693	0
        CPPUNIT_ASSERT_EQUAL ( (int)121693, this->txtReader_->position_[1][9] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->info_[23-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->content_[23-2][0] );

        //Pf3D7_04_v3	149598	0
        CPPUNIT_ASSERT_EQUAL ( (int)149598, this->txtReader_->position_[3][4] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->info_[47-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->content_[47-2][0] );

        //Pf3D7_06_v3	180217	0
        CPPUNIT_ASSERT_EQUAL ( (int)180217, this->txtReader_->position_[5][5] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->info_[91-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->content_[91-2][0] );

        //Pf3D7_06_v3	180270	0
        CPPUNIT_ASSERT_EQUAL ( (int)180270, this->txtReader_->position_[5][8] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->info_[94-2] );
        CPPUNIT_ASSERT_EQUAL ( (double)0, this->txtReader_->content_[94-2][0] );

        this->afterExclude_->findAndKeepMarkers (excludedMarkers_);
        CPPUNIT_ASSERT_EQUAL ( (size_t)93, this->afterExclude_->info_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)93, this->afterExclude_->nLoci_ );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->afterExclude_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)6, this->afterExclude_->position_.size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)12, this->afterExclude_->position_[0].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)(16-4), this->afterExclude_->position_[1].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)17, this->afterExclude_->position_[2].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)(24-1), this->afterExclude_->position_[3].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)20, this->afterExclude_->position_[4].size() );
        CPPUNIT_ASSERT_EQUAL ( (size_t)(11-2), this->afterExclude_->position_[5].size() );

        CPPUNIT_ASSERT_EQUAL (this->txtReader_->info_.size(), this->afterExclude_->info_.size() );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->content_.size(), this->afterExclude_->content_.size() );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->keptContent_.size(), this->afterExclude_->keptContent_.size() );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->keptContent_.size(), (size_t)0 );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->nInfoLines_, this->afterExclude_->nInfoLines_ );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->nInfoLines_, (size_t)1);

        for ( size_t i = 0; i < 93; i++) {
            CPPUNIT_ASSERT_EQUAL (this->txtReader_->info_[i], this->afterExclude_->info_[i] );
            CPPUNIT_ASSERT_EQUAL (this->txtReader_->content_[i][0], this->afterExclude_->content_[i][0] );
        }


        for ( size_t i = 0; i < 93; i++) {
            CPPUNIT_ASSERT_EQUAL (this->txtReader_->info_[i], this->afterExclude_->info_[i] );
            CPPUNIT_ASSERT_EQUAL (this->txtReader_->content_[i][0], this->afterExclude_->content_[i][0] );
        }

        CPPUNIT_ASSERT_EQUAL (this->txtReader_->chrom_.size(), this->afterExclude_->chrom_.size() );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->chrom_.size(), (size_t)6 );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->indexOfChromStarts_.size(), this->afterExclude_->indexOfChromStarts_.size() );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->indexOfChromStarts_.size(), (size_t)6 );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->position_.size(), this->afterExclude_->position_.size() );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->position_.size(), (size_t)6 );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->keptPosition_.size(), this->afterExclude_->keptPosition_.size() );
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->keptPosition_.size(), (size_t)0 );

        for ( size_t i = 0; i < 6; i++) {
            CPPUNIT_ASSERT_EQUAL (this->txtReader_->chrom_[i], this->afterExclude_->chrom_[i] );
            CPPUNIT_ASSERT_EQUAL (this->txtReader_->indexOfChromStarts_[i], this->afterExclude_->indexOfChromStarts_[i] );
            CPPUNIT_ASSERT_EQUAL (this->txtReader_->position_[i].size(), this->afterExclude_->position_[i].size() );
            for ( size_t j = 0; j < this->txtReader_->position_[i].size(); j++) {
                CPPUNIT_ASSERT_EQUAL (this->txtReader_->position_[i][j], this->afterExclude_->position_[i][j] );

            }
        }
        CPPUNIT_ASSERT_EQUAL (this->txtReader_->nLoci_, this->afterExclude_->nLoci_ );
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestTxtReader );
