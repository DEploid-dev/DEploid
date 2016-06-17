#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "pfDeconvIO.hpp"

class TestIO : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestIO );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testInitialization );
    CPPUNIT_TEST( testPrintHelp );
    CPPUNIT_TEST( testNotEnoughArg );
    CPPUNIT_TEST( testWrongType );
    CPPUNIT_TEST( testFileNameMissing );
    CPPUNIT_TEST( testInvalidInputFile );
    CPPUNIT_TEST( testUnknowArg );
    CPPUNIT_TEST( testFlagsConflict );
    CPPUNIT_TEST( testExtractRefAltPlaf );
    CPPUNIT_TEST( testLociNumberUnequal );
    CPPUNIT_TEST( testForbidMoves );
    CPPUNIT_TEST( testInitialProp );
    CPPUNIT_TEST( testExtractPostProbException );
    CPPUNIT_TEST_SUITE_END();

  private:
    PfDeconvIO* input_;
    double epsilon3;

  public:
    void setUp() {
        this->epsilon3 = 0.000000000001;
        this->input_ = new PfDeconvIO ();
    }

    void tearDown() {
        delete input_;
    }

    void testInitialization(){
        CPPUNIT_ASSERT_EQUAL( this->input_->randomSeedWasSet_, false );
        CPPUNIT_ASSERT_EQUAL( this->input_->initialPropWasGiven(), false );
        CPPUNIT_ASSERT_EQUAL( this->input_->exclude_sites_ , false );
        CPPUNIT_ASSERT( this->input_->excludedMarkers == NULL);
        CPPUNIT_ASSERT_EQUAL( this->input_->random_seed_, (size_t)0 );
        CPPUNIT_ASSERT_EQUAL( this->input_->help(), false);
        CPPUNIT_ASSERT_EQUAL( this->input_->usePanel(), true);
        CPPUNIT_ASSERT_EQUAL( this->input_->precision_ , (size_t)8);
        CPPUNIT_ASSERT( this->input_->prefix_ == "pf3k-pfDeconv");
        CPPUNIT_ASSERT_EQUAL( this->input_->kStrain_ , (size_t)5);
        CPPUNIT_ASSERT_EQUAL( this->input_->nMcmcSample_ , (size_t)800);
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdateProp() , true);
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdatePair() , true);
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdateSingle() , true);
        CPPUNIT_ASSERT_EQUAL( this->input_->doExportPostProb() , false);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( this->input_->mcmcBurn_, 0.5, epsilon3);
        CPPUNIT_ASSERT_EQUAL( this->input_->mcmcMachineryRate_ , (size_t)5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( this->input_->missCopyProb_ , 0.01, epsilon3);
        CPPUNIT_ASSERT_EQUAL( this->input_->useConstRecomb() , false );
        CPPUNIT_ASSERT_EQUAL( this->input_->forbidCopyFromSame() , false );
        CPPUNIT_ASSERT_DOUBLES_EQUAL( this->input_->constRecombProb_ , 1.0, epsilon3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL( this->input_->averageCentimorganDistance_, 15000.0, epsilon3 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL( this->input_->Ne_ , 10.0, epsilon3 );
    }

    void testMainConstructor() {

        char *argv[] = { "./pfDeconv" };
        CPPUNIT_ASSERT_NO_THROW ( this->input_->core(1, argv) );

        CPPUNIT_ASSERT_EQUAL((size_t)0, input_->argv_.size());
        CPPUNIT_ASSERT_EQUAL( true, input_->help());

        char *argv1[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_NO_THROW ( this->input_->core(11, argv1) );
    }


    void testInitialProp() {

        char *argv1[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-initialP", "0.1", "0.2", "0.3", "0.4", "-k", "4" };
        CPPUNIT_ASSERT_NO_THROW ( this->input_->core(16, argv1) );
        CPPUNIT_ASSERT_THROW ( this->input_->core(14, argv1), NumOfPropNotMatchNumStrain );

        char *argv2[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-initialP", "0.1", "0.2", "0.3", "0.3", "0.2"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(15, argv2), SumOfPropNotOne );

        char *argv3[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-initialP", "0.1", "0.2", "0.3", "0.3", "0.02"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(15, argv3), SumOfPropNotOne );

        char *argv4[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-initialP", "-o", "tmp"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(12, argv4), NotEnoughArg );

        char *argv5[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-initialP"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv5), NotEnoughArg );
    }


    void testPrintHelp(){
        CPPUNIT_ASSERT_NO_THROW ( this->input_->printHelp() );
    }


    void testNotEnoughArg(){
        char *argv1[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-o" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv1), NotEnoughArg );

        char *argv2[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-seed"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv2), NotEnoughArg );

        char *argv3[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-p"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv3), NotEnoughArg );

        char *argv4[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-k"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv4), NotEnoughArg );

        char *argv5[] = { "./pfDeconv",
                         "-ref",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(8, argv5), NotEnoughArg );

        char *argv6[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(8, argv6), NotEnoughArg );

        char *argv7[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf",
                         "-panel", "tests/testData/lab_first100_Panel.txt" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(8, argv7), NotEnoughArg );

        char *argv8[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(8, argv8), NotEnoughArg );
    }


    void testWrongType(){
        char *argv2[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-seed", "asdf"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv2), NotEnoughArg );
        CPPUNIT_ASSERT_THROW ( this->input_->core(11, argv2), WrongType );

        char *argv3[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-p", "asdf"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv3), NotEnoughArg );
        CPPUNIT_ASSERT_THROW ( this->input_->core(11, argv3), WrongType );

        char *argv4[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-k", "asdf"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv4), NotEnoughArg );
        CPPUNIT_ASSERT_THROW ( this->input_->core(11, argv4), WrongType );

        char *argv5[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-nSample", "asdf"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv5), NotEnoughArg );
        CPPUNIT_ASSERT_THROW ( this->input_->core(11, argv5), WrongType );

        char *argv6[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-rate", "asdf"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv6), NotEnoughArg );
        CPPUNIT_ASSERT_THROW ( this->input_->core(11, argv6), WrongType );
    }


    void testFileNameMissing(){
        char *argv1[] = { "./pfDeconv",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(9, argv1), FileNameMissing );
        CPPUNIT_ASSERT_THROW_MESSAGE( "Ref count file path missing!", this->input_->core(9, argv1), FileNameMissing );
        try {
            this->input_->core(9, argv1);
        }
        catch (const exception &e) {
            CPPUNIT_ASSERT_EQUAL( string("Ref count file path missing!"), string(e.what()) );
        }

        char *argv2[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(9, argv2), FileNameMissing );

        try {
            this->input_->core(9, argv2);
        }
        catch (const exception &e) {
            CPPUNIT_ASSERT_EQUAL( string("Alt count file path missing!"), string(e.what()) );
        }

        char *argv3[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(9, argv3), FileNameMissing );

        char *argv4[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(9, argv4), FileNameMissing );
    }


    void testInvalidInputFile(){
        char *argv[] = { "./pfDeconv",
                         "-ref", "PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-o", "tmp1" };
        CPPUNIT_ASSERT_THROW( this->input_->core(11, argv), InvalidInputFile );
        try {
            this->input_->core(11, argv);
        }
        catch (const exception &e) {
            CPPUNIT_ASSERT_EQUAL( string("Invalid input file: PG0390_first100ref.txt"), string(e.what()) );
        }
    }


    void testUnknowArg(){
        char *argv[] = { "./pfDeconv", "-unknow"};
        CPPUNIT_ASSERT_THROW( this->input_->core(2, argv), UnknowArg );
    }


    void testFlagsConflict(){
        char *argv1[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-noPanel" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv1), FlagsConflict );

        char *argv2[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-noPanel",
                         "-panel", "tests/testData/lab_first100_Panel.txt"};
        CPPUNIT_ASSERT_THROW ( this->input_->core(10, argv2), FlagsConflict );
    }


    void testLociNumberUnequal(){
        char *argv1[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_samples_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(9, argv1), LociNumberUnequal );

        char *argv2[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390.C_alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt" };
        CPPUNIT_ASSERT_THROW ( this->input_->core(9, argv2), LociNumberUnequal );

    }


    void testExtractRefAltPlaf(){
        char *argv1[] = { "./pfDeconv",
         "-ref", "tests/testData/PG0390_first100ref.txt",
         "-alt", "tests/testData/PG0390_first100alt.txt",
         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
         "-panel", "tests/testData/lab_first100_Panel.txt",
         "-o", "tmp1" };
        CPPUNIT_ASSERT_NO_THROW ( this->input_->core(11, argv1) );

        CPPUNIT_ASSERT_EQUAL((size_t)100, this->input_->nLoci_);
        //"Pf3D7_01_v3"	93157	0.0190612159917058
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0190612159917058, this->input_->plaf_[0], 0.000000000001);
        //"Pf3D7_01_v3"	93157	85
        CPPUNIT_ASSERT_EQUAL(85.0, this->input_->refCount_[0]);
        //"Pf3D7_01_v3"	93157	0
        CPPUNIT_ASSERT_EQUAL(0.0, this->input_->altCount_[0]);

        //"Pf3D7_01_v3"	95518	0.687463394087723
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.687463394087723, this->input_->plaf_[4], 0.000000000001);
        //"Pf3D7_01_v3"	95518	156
        CPPUNIT_ASSERT_EQUAL(156.0, this->input_->refCount_[4]);
        //"Pf3D7_01_v3"	95518	46
        CPPUNIT_ASSERT_EQUAL(46.0, this->input_->altCount_[4]);

        //"Pf3D7_01_v3"	113396	0.0207016179419847
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0207016179419847, this->input_->plaf_[20], 0.000000000001);
        //"Pf3D7_01_v3"	113396	177
        CPPUNIT_ASSERT_EQUAL(177.0, this->input_->refCount_[20]);
        //"Pf3D7_01_v3"	113396	0
        CPPUNIT_ASSERT_EQUAL(0.0, this->input_->altCount_[20]);

        //"Pf3D7_01_v3"	180270	0.426732645350475
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.426732645350475, this->input_->plaf_[99], 0.000000000001);
        //"Pf3D7_01_v3"	180270	199
        CPPUNIT_ASSERT_EQUAL(199.0, this->input_->refCount_[99]);
        //"Pf3D7_01_v3"	180270	0
        CPPUNIT_ASSERT_EQUAL(0.0, this->input_->altCount_[99]);
    }


    void testForbidMoves() {
        char *argv1[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-o", "tmp1", "-forbidUpdateProp", "-forbidUpdateSingle", "-forbidUpdatePair" };
        CPPUNIT_ASSERT_NO_THROW ( this->input_->core(14, argv1) );
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdateProp(), false);
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdateSingle(), false);
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdatePair(), false);

        CPPUNIT_ASSERT_NO_THROW ( this->input_->core(13, argv1) );
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdateProp(), false);
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdateSingle(), false);
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdatePair(), true);

        CPPUNIT_ASSERT_NO_THROW ( this->input_->core(12, argv1) );
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdateProp(), false);
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdateSingle(), true);
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdatePair(), true);

        CPPUNIT_ASSERT_NO_THROW ( this->input_->core(11, argv1) );
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdateProp() , true);
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdateSingle() , true);
        CPPUNIT_ASSERT_EQUAL( this->input_->doUpdatePair() , true);
    }

    void testExtractPostProbException(){
        char *argv1[] = { "./pfDeconv",
                         "-ref", "tests/testData/PG0390_first100ref.txt",
                         "-alt", "tests/testData/PG0390_first100alt.txt",
                         "-plaf", "tests/testData/labStrains_first100_PLAF.txt",
                         "-panel", "tests/testData/lab_first100_Panel.txt",
                         "-o", "tmp1", "-exportPostProb", "-k", "2" };
        CPPUNIT_ASSERT_NO_THROW ( this->input_->core(14, argv1) );
        CPPUNIT_ASSERT_EQUAL( this->input_->doExportPostProb() , true);
        CPPUNIT_ASSERT_THROW ( this->input_->core(12, argv1), onlyExportPostProbWhenTwoStrains );

    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestIO );
