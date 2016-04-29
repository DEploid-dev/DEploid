#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "updateHap.hpp"

class TestUpdatePairHap : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestUpdatePairHap );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testComputeMarginalDist );
    CPPUNIT_TEST_SUITE_END();

  private:
    UpdatePairHap * updatePairHapPanel_;
    UpdatePairHap * updatePairHapPlaf_;
    vector <double> refCount_;
    vector <double> altCount_;
    vector <double> plaf_;
    vector <double> expectedWsaf_;
    vector <double> proportion_;
    vector < vector <double> > haplotypes_;
    MersenneTwister* rg_;
    size_t segmentStartIndex_;
    size_t nLoci_;
    Panel* panel_;
    double missCopyProb_;
    size_t strainIndex1_;
    size_t strainIndex2_;
    bool forbidCopyFromSame_ ;

  public:
    void setUp(){
        this->expectedWsaf_ = vector <double> ();
        this->rg_ = new MersenneTwister((size_t)1);
        this->panel_ = new Panel();
        this->nLoci_ = 0;
        this->forbidCopyFromSame_ = true;


        this->updatePairHapPlaf_ = new UpdatePairHap( refCount_,
                                          altCount_,
                                          plaf_,
                                          expectedWsaf_,
                                          proportion_,
                                          haplotypes_,
                                          rg_,
                                          segmentStartIndex_,
                                          nLoci_,
                                          NULL,
                                          missCopyProb_, forbidCopyFromSame_,
                                          strainIndex1_, strainIndex2_ );

        this->updatePairHapPanel_ = new UpdatePairHap( refCount_,
                                          altCount_,
                                          plaf_,
                                          expectedWsaf_,
                                          proportion_,
                                          haplotypes_,
                                          rg_,
                                          segmentStartIndex_,
                                          nLoci_,
                                          panel_,
                                          missCopyProb_, forbidCopyFromSame_,
                                          strainIndex1_, strainIndex2_ );
    }


    void tearDown(){
        delete updatePairHapPlaf_;
        delete updatePairHapPanel_;
        delete panel_;
        this->rg_->clearFastFunc();
        delete rg_;
    }

    void testMainConstructor(){ }

    void testComputeMarginalDist(){
        vector < vector < double> > tmpMat;
        tmpMat.push_back( vector <double> ({1.0, 2.0, 3.0}) );
        tmpMat.push_back( vector <double> ({4.0, 5.0, 6.0}) );
        tmpMat.push_back( vector <double> ({7.0, 8.0, 9.0}) );

        vector <double> rowMarg = this->updatePairHapPanel_->computeRowMarginalDist(tmpMat);
        CPPUNIT_ASSERT_EQUAL((size_t)3, rowMarg.size());
        CPPUNIT_ASSERT_EQUAL(6.0,  rowMarg[0]);
        CPPUNIT_ASSERT_EQUAL(15.0, rowMarg[1]);
        CPPUNIT_ASSERT_EQUAL(24.0, rowMarg[2]);

        vector <double> colMarg = this->updatePairHapPanel_->computeColMarginalDist(tmpMat);
        CPPUNIT_ASSERT_EQUAL((size_t)3, colMarg.size());
        CPPUNIT_ASSERT_EQUAL(12.0, colMarg[0]);
        CPPUNIT_ASSERT_EQUAL(15.0, colMarg[1]);
        CPPUNIT_ASSERT_EQUAL(18.0, colMarg[2]);
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestUpdatePairHap );
