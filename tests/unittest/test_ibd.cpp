#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "src/ibd.hpp"
#include <iomanip>      // std::setw


class TestIBDUtility: public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestIBDUtility );
    CPPUNIT_TEST( testConvertIntToBinary );
    CPPUNIT_TEST( testNchoose2 );
    CPPUNIT_TEST_SUITE_END();

    void testConvertIntToBinary(){
        CPPUNIT_ASSERT_THROW (convertIntToBinary(20, 3), OutOfVectorSize);
        CPPUNIT_ASSERT_THROW (convertIntToBinary(15, 3), OutOfVectorSize);
        CPPUNIT_ASSERT_THROW (convertIntToBinary(8, 3), OutOfVectorSize);
        CPPUNIT_ASSERT_THROW_MESSAGE( "Out of vector size!", convertIntToBinary(8, 3), OutOfVectorSize);

        CPPUNIT_ASSERT_NO_THROW (convertIntToBinary(7, 3));

        vector <int> test1 = convertIntToBinary(3, 4);
        CPPUNIT_ASSERT_EQUAL (test1.size(), (size_t)4);
        CPPUNIT_ASSERT_EQUAL (test1[0], (int)0);
        CPPUNIT_ASSERT_EQUAL (test1[1], (int)0);
        CPPUNIT_ASSERT_EQUAL (test1[2], (int)1);
        CPPUNIT_ASSERT_EQUAL (test1[3], (int)1);

        CPPUNIT_ASSERT_NO_THROW(enumerateBinaryMatrixOfK(10));
        vector < vector <int> > test2 = enumerateBinaryMatrixOfK(4);
        CPPUNIT_ASSERT_EQUAL (test2.size(), (size_t)16);
        for (size_t i = 0; i < 16; i++){
            CPPUNIT_ASSERT_EQUAL (test2[i].size(), (size_t)4);
            CPPUNIT_ASSERT_EQUAL ((int)(8*test2[i][0] +
                                        4*test2[i][1] +
                                        2*test2[i][2] +
                                        1*test2[i][3]), (int)i);
        }
    }

    void testNchoose2(){
        CPPUNIT_ASSERT_THROW (nchoose2(-1), InvalidInput);
        CPPUNIT_ASSERT_THROW (nchoose2(0), InvalidInput);
        CPPUNIT_ASSERT_THROW (nchoose2(1), InvalidInput);
        CPPUNIT_ASSERT_NO_THROW (nchoose2(2));
        CPPUNIT_ASSERT_EQUAL (nchoose2(2), (int)1);
        CPPUNIT_ASSERT_EQUAL (nchoose2(3), (int)3);
        CPPUNIT_ASSERT_EQUAL (nchoose2(4), (int)6);
        CPPUNIT_ASSERT_EQUAL (nchoose2(5), (int)10);
    }
};


class TestIBDconfig : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestIBDconfig );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testCountEffectiveK );
    CPPUNIT_TEST_SUITE_END();

  private:
    IBDconfiguration* ibd3_;
    IBDconfiguration* ibd4_;
    IBDconfiguration* ibd5_;

  public:
    void setUp() {
        ibd3_ = new IBDconfiguration();
        CPPUNIT_ASSERT_NO_THROW( ibd3_->buildIBDconfiguration(3));
        ibd4_ = new IBDconfiguration();
        CPPUNIT_ASSERT_NO_THROW( ibd4_->buildIBDconfiguration(4));
        ibd5_ = new IBDconfiguration();
        CPPUNIT_ASSERT_NO_THROW( ibd5_->buildIBDconfiguration(5));
    }


    void tearDown() {
        delete ibd3_;
        delete ibd4_;
        delete ibd5_;
    }


    vector <int> computeEffectiveKCount(int k, IBDconfiguration* ibd){
        vector <int> ret(k);
        for ( int index = 0; index < k; index++){
            for (int effk : ibd->effectiveK ){
                if (effk == (index+1)){
                    ret[index]++;
                }
            }
        }
        return ret;
    }



    void testMainConstructor(){
        CPPUNIT_ASSERT_NO_THROW( IBDconfiguration());
        //CPPUNIT_ASSERT_NO_THROW( IBDconfiguration(5));
        CPPUNIT_ASSERT_EQUAL( ibd3_->states.size(), (size_t)5);
        CPPUNIT_ASSERT_EQUAL( ibd4_->states.size(), (size_t)15);
        CPPUNIT_ASSERT_EQUAL( ibd5_->states.size(), (size_t)52);
    }


    void testCountEffectiveK(){
        // k = 3
        vector <int> count3 = computeEffectiveKCount(3, ibd3_);
        //table(make.ibd.mat.joe(3)$k.eff)
        //1 2 3
        //1 3 1
        CPPUNIT_ASSERT_EQUAL( count3[0], (int)1); // 1
        CPPUNIT_ASSERT_EQUAL( count3[1], (int)3); // 2
        CPPUNIT_ASSERT_EQUAL( count3[2], (int)1); // 3

        // k = 4
        vector <int> count4 = computeEffectiveKCount(4, ibd4_);
        //table(make.ibd.mat.joe(4)$k.eff)
        //1 2 3 4
        //1 7 6 1
        CPPUNIT_ASSERT_EQUAL( count4[0], (int)1); // 1
        CPPUNIT_ASSERT_EQUAL( count4[1], (int)7); // 2
        CPPUNIT_ASSERT_EQUAL( count4[2], (int)6); // 3
        CPPUNIT_ASSERT_EQUAL( count4[3], (int)1); // 4

        // k = 5
        vector <int> count5 = computeEffectiveKCount(5, ibd5_);
        //table(make.ibd.mat.joe(5)$k.eff)
        //1  2  3  4  5
        //1 15 25 10  1
        CPPUNIT_ASSERT_EQUAL( count5[0], (int)1);  // 1
        CPPUNIT_ASSERT_EQUAL( count5[1], (int)15); // 2
        CPPUNIT_ASSERT_EQUAL( count5[2], (int)25); // 3
        CPPUNIT_ASSERT_EQUAL( count5[3], (int)10); // 4
        CPPUNIT_ASSERT_EQUAL( count5[4], (int)1);  // 5
    }


};

class TestHprior : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestHprior );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testCheckColSums);
    CPPUNIT_TEST( testCheckHSetColSums);
    CPPUNIT_TEST_SUITE_END();

  private:
    Hprior* hprior3_;
    Hprior* hprior4_;
    Hprior* hprior5_;
    vector <double> plaf;
  public:
    void setUp() {
        plaf = vector <double> ({0.0, 0.5, 1.0});
        hprior3_ = new Hprior;
        CPPUNIT_ASSERT_NO_THROW(hprior3_->buildHprior(3, plaf));
        hprior4_ = new Hprior;(4, plaf);
        CPPUNIT_ASSERT_NO_THROW(hprior4_->buildHprior(4, plaf));
        hprior5_ = new Hprior;(5, plaf);
        CPPUNIT_ASSERT_NO_THROW(hprior5_->buildHprior(5, plaf));
    }


    void tearDown() {
        delete hprior3_;
        delete hprior4_;
        delete hprior5_;
    }


    void testMainConstructor(){
        CPPUNIT_ASSERT_EQUAL( hprior3_->nState(), (size_t)22); // This is equal to sum(2^IBDconfiguration(3).effectiveK)
        CPPUNIT_ASSERT_EQUAL( hprior4_->nState(), (size_t)94);
        CPPUNIT_ASSERT_EQUAL( hprior5_->nState(), (size_t)454);
    }


    void testCheckColSumsCore(Hprior* hprior){
            //cout <<endl;
        //for (vector <double> prob : hprior->priorProb){
            //for ( double p : prob ){
                 //cout << setw(10) << p;
            //}
            //cout <<endl;
            ////cout << "i = "<<i<<" tmpSum = "<< tmpSum<<endl;
            ////CPPUNIT_ASSERT_EQUAL(tmpSum, (double)hprior->nPattern());
        //}

        for ( size_t i = 0; i < hprior->nLoci(); i++){
            double tmpSum = 0.0;
            for (vector <double> prob : hprior->priorProb){
                tmpSum += prob[i];
            }
            //cout << "i = "<<i<<" tmpSum = "<< tmpSum<<endl;
            CPPUNIT_ASSERT_EQUAL(tmpSum, (double)hprior->nPattern());
        }
    }


    void testCheckColSums(){
        this->testCheckColSumsCore(hprior3_);
        this->testCheckColSumsCore(hprior4_);
        this->testCheckColSumsCore(hprior5_);
    }


    void computeHSetColSumsCore(Hprior *hprior, int sum){
        for ( size_t i = 0; i < hprior->kStrain(); i++){
            int tmpSum = 0;
            for (vector <int> hRow : hprior->hSet){
                tmpSum += hRow[i];
            }
            CPPUNIT_ASSERT_EQUAL(tmpSum, sum);
        }
    }


    void testCheckHSetColSums(){
        this->computeHSetColSumsCore(hprior3_, 11);
        this->computeHSetColSumsCore(hprior4_, 47);
        this->computeHSetColSumsCore(hprior5_, 227);
    }
};


class TestIBDpath : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestIBDpath );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testmakeLlkSurf );
    CPPUNIT_TEST( testIbdTransProbs );
    CPPUNIT_TEST( testComputeUniqueEffectiveKCount );
    CPPUNIT_TEST( testIBDconfigureHeader );
    CPPUNIT_TEST( testStatePrior );
    CPPUNIT_TEST( checkFwd );
    CPPUNIT_TEST( checkBwd );
    CPPUNIT_TEST( checkPainting );
    CPPUNIT_TEST_SUITE_END();

  private:
    IBDpath* ibdPath2_;
    IBDpath* ibdPath3_;
    IBDpath* ibdPath5_;
    DEploidIO* dEploidIO_;
    MersenneTwister* rg_;
    double epsilon2;

  public:
    void setUp() {
        dEploidIO_ = new DEploidIO();
        dEploidIO_->altCount_ = vector<double> ({2, 100, 50, 50, 2, 2});
        dEploidIO_->refCount_ = vector<double> ({100, 2, 0, 50, 0, 2});
        dEploidIO_->setKstrain(3);
        dEploidIO_->plaf_ = vector<double> ({0.1, .4, .4, .3, .2, .5});
        vector < vector <int> > testPosition;
        testPosition.push_back(vector<int> ({1,2,3,4,5,6}));
        //testPosition.push_back(vector<int> ({200, 3000}));
        //testPosition.push_back(vector<int> ({300}));
        //testPosition.push_back(vector<int> ({300, 400, 500}));
        dEploidIO_->position_ = testPosition;
        dEploidIO_->nLoci_ = 6;
        //dEploidIO_->chrom_ = vector <string> ({"chrom1", "chrom2", "chrom3"});
        dEploidIO_->chrom_ = vector <string> ({"chrom1"});
        dEploidIO_->useConstRecomb_ = true;
        dEploidIO_->constRecombProb_ = 0.000001;
        rg_ = new MersenneTwister(dEploidIO_->randomSeed());
        epsilon2 = 0.001;
        ibdPath3_ = new IBDpath;
        this->ibdPath3_->init(*dEploidIO_, rg_);
        dEploidIO_->setKstrain(5);
        ibdPath5_ = new IBDpath;
        this->ibdPath5_->init(*dEploidIO_, rg_);
        dEploidIO_->setKstrain(2);
        ibdPath2_ = new IBDpath;
        this->ibdPath2_->init(*dEploidIO_, rg_);
    }


    void tearDown() {
        delete ibdPath2_;
        delete ibdPath3_;
        delete ibdPath5_;
        delete rg_;
        delete dEploidIO_;
    }

    void testMainConstructor(){
        //for (double p : ibdPath2_->ibdRecombProbs.pRec_){cout <<p<<endl;}
        //for (double p : ibdPath3_->ibdRecombProbs.pRec_){cout <<p<<endl;}
    }

    void testmakeLlkSurf(){
        //> make.llk.surf(c(2,100,50,50,2,2), c(102,102,50,100,2,4), do.plot=FALSE)
                  //[,1]      [,2]
        //[1,]  2.149687 62.602446
        //[2,] 62.602446  2.149687
        //[3,] 50.250071  1.507423
        //[4,] 24.735000 24.735000
        //[5,]  2.968471  1.029853
        //[6,]  2.800135  2.800135

        CPPUNIT_ASSERT_EQUAL( this->ibdPath3_->llkSurf.size(), (size_t)6 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.149687, this->ibdPath3_->llkSurf[0][0], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 62.602446, this->ibdPath3_->llkSurf[0][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 62.602446, this->ibdPath3_->llkSurf[1][0], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.149687, this->ibdPath3_->llkSurf[1][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 50.250071, this->ibdPath3_->llkSurf[2][0], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1.507423, this->ibdPath3_->llkSurf[2][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 24.735000, this->ibdPath3_->llkSurf[3][0], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 24.735000, this->ibdPath3_->llkSurf[3][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.968471, this->ibdPath3_->llkSurf[4][0], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1.029853, this->ibdPath3_->llkSurf[4][1], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.800135, this->ibdPath3_->llkSurf[5][0], epsilon2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.800135, this->ibdPath3_->llkSurf[5][1], epsilon2);
    }


    void testIbdTransProbs(){
        //k = 3
        //> tij
             //[,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
        //[1,]    1    1    1    1    1    1    1    1    0     0     0     0     0     0
        //[2,]    0    0    0    0    0    0    0    0    1     1     1     1     0     0
        //[3,]    0    0    0    0    0    0    0    0    0     0     0     0     1     1
        //[4,]    0    0    0    0    0    0    0    0    0     0     0     0     0     0
        //[5,]    0    0    0    0    0    0    0    0    0     0     0     0     0     0
             //[,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
        //[1,]     0     0     0     0     0     0     0     0
        //[2,]     0     0     0     0     0     0     0     0
        //[3,]     1     1     0     0     0     0     0     0
        //[4,]     0     0     1     1     1     1     0     0
        //[5,]     0     0     0     0     0     0     1     1
        //vector <double> tmpPlaf = vector <double> ({0.1, 0.2, 0.3});
        //CPPUNIT_ASSERT_NO_THROW(this->ibdPath3_->hprior.buildHprior(3, tmpPlaf));
        CPPUNIT_ASSERT_EQUAL(this->ibdPath3_->hprior.nPattern(), (size_t)5);
        CPPUNIT_ASSERT_EQUAL(this->ibdPath3_->hprior.nState(), (size_t)22);

        //CPPUNIT_ASSERT_NO_THROW(this->ibdPath3_->makeIbdTransProbs());
        CPPUNIT_ASSERT_EQUAL((size_t)5, this->ibdPath3_->ibdTransProbs.size());
        double tmpValue = sumOfMat(this->ibdPath3_->ibdTransProbs);
        CPPUNIT_ASSERT_EQUAL((double)22.0, tmpValue);
        for ( size_t i =  0; i < 22; i++){
            CPPUNIT_ASSERT_EQUAL(this->ibdPath3_->ibdTransProbs[this->ibdPath3_->hprior.stateIdx[i]][i], (double)1.0);
        }
    }


    void testComputeUniqueEffectiveKCount(){
        //vector <double> tmpPlaf = vector <double> ({0.1, 0.2, 0.3});
        //CPPUNIT_ASSERT_NO_THROW(this->ibdPath3_->hprior.buildHprior(5, tmpPlaf));
        CPPUNIT_ASSERT_NO_THROW(this->ibdPath5_->computeUniqueEffectiveKCount());
        CPPUNIT_ASSERT_EQUAL(this->ibdPath5_->uniqueEffectiveKCount.size(), (size_t)5);
        //> table(make.ibd.mat.joe(5)$k.eff)
        //1  2  3  4  5
        //1 15 25 10  1
        CPPUNIT_ASSERT_EQUAL(this->ibdPath5_->uniqueEffectiveKCount[0], (int)1);
        CPPUNIT_ASSERT_EQUAL(this->ibdPath5_->uniqueEffectiveKCount[1], (int)15);
        CPPUNIT_ASSERT_EQUAL(this->ibdPath5_->uniqueEffectiveKCount[2], (int)25);
        CPPUNIT_ASSERT_EQUAL(this->ibdPath5_->uniqueEffectiveKCount[3], (int)10);
        CPPUNIT_ASSERT_EQUAL(this->ibdPath5_->uniqueEffectiveKCount[4], (int)1);
    }


    void testIBDconfigureHeader(){
        vector <string> headerOf3Strain = this->ibdPath3_->getIBDprobsHeader();
        CPPUNIT_ASSERT_EQUAL((size_t)5, headerOf3Strain.size());
        CPPUNIT_ASSERT("0-1-2" == headerOf3Strain[0]);
        CPPUNIT_ASSERT("0-2-2" == headerOf3Strain[1]);
        CPPUNIT_ASSERT("2-1-2" == headerOf3Strain[2]);
        CPPUNIT_ASSERT("2-2-2" == headerOf3Strain[3]);
        CPPUNIT_ASSERT("1-1-2" == headerOf3Strain[4]);
    }


    void testStatePrior(){
        //vector <double> effectiveKPrior = this->ibdPath3_->computeEffectiveKPrior(this->ibdPath3_->theta());
        vector <double> effectiveKPrior = vector <double> (this->ibdPath3_->hprior.nPattern(), 1.0/this->ibdPath3_->hprior.nPattern());
        //for (double p : effectiveKPrior){cout<<p<<endl;}
        vector <double> statePrior = this->ibdPath3_->computeStatePrior(effectiveKPrior);
        //for (double p : statePrior){cout<<p<<endl;}
    }

    void checkFwd(){
        vector <double> effectiveKPrior = vector <double> (this->ibdPath2_->hprior.nPattern(), 1.0/this->ibdPath2_->hprior.nPattern());
        //for (double p : effectiveKPrior){cout<<p<<endl;}
        vector <double> statePrior = this->ibdPath2_->computeStatePrior(effectiveKPrior);
        this->ibdPath2_->computeIbdPathFwdProb(vector <double> ({.5, .5}), statePrior);
        vector <double > tmp =this->ibdPath2_->fm[5];
        normalizeBySum(tmp);
        //for (double p:tmp) printf("%8.2f\n", p);
        vector < vector <double>> reshapedFwd = this->ibdPath2_->reshapeProbs(this->ibdPath2_->fm);
        for ( vector <double> p_vec : reshapedFwd){
            for (double p : p_vec){
                cout<<p<<" ";
            }
            cout<<endl;
        }


    }

    void checkBwd(){
        vector <double> effectiveKPrior = vector <double> (this->ibdPath2_->hprior.nPattern(), 1.0/this->ibdPath2_->hprior.nPattern());
        //for (double p : effectiveKPrior){cout<<p<<endl;}
        vector <double> statePrior = this->ibdPath2_->computeStatePrior(effectiveKPrior);
        this->ibdPath2_->computeIbdPathBwdProb(vector <double> ({.5, .5}), effectiveKPrior, statePrior);

        for (size_t i = 0; i < 6; i++){
            for (size_t j = 0; j < 6; j++){
                printf("%8.3f  ", this->ibdPath2_->bwd[j][i]);
            }
            cout << endl;
        }

            vector < vector <double>> reshapedBwd = this->ibdPath2_->reshapeProbs(this->ibdPath2_->bwd);
        for ( vector <double> p_vec : reshapedBwd){
            for (double p : p_vec){
                cout<<p<<" ";
            }
            cout<<endl;
        }

        //for (size_t i = 0; i < 6; i++){
            //vector <double > tmp = this->ibdPath2_->bwd[i];
            //normalizeBySum(tmp);
        //}
    }

    void checkPainting(){
        CPPUNIT_ASSERT_NO_THROW(this->ibdPath2_->buildPathProbabilityForPainting(vector <double> ({.5, .5})));
        CPPUNIT_ASSERT_NO_THROW(this->ibdPath3_->buildPathProbabilityForPainting(vector <double> ({.5, .25, .25})));
        CPPUNIT_ASSERT_NO_THROW(this->ibdPath5_->buildPathProbabilityForPainting(vector <double> ({.25, .2, .2, .2, .15})));
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION(TestIBDUtility);
CPPUNIT_TEST_SUITE_REGISTRATION(TestIBDconfig);
CPPUNIT_TEST_SUITE_REGISTRATION(TestHprior);
CPPUNIT_TEST_SUITE_REGISTRATION(TestIBDpath);
