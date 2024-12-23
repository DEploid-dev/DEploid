/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of dEploid.
 *
 * dEploid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <ctime>
#include <iterator>
#include <cassert>        // assert
#include "utility.hpp"    // normailize by sum
#include "dEploidIO.hpp"
#include "ibd.hpp"
#include "lasso/src/dEploidLasso.hpp"

DEploidIO::DEploidIO() {
    this->init();
}


DEploidIO::DEploidIO(const std::string &arg) {
    this->init();
    std::istringstream iss(arg);
    copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter(argv_));
    this->argv_i = argv_.begin();
    this->core();
}


DEploidIO::DEploidIO(int argc, char *argv[]) {
    this->init();
    argv_ = std::vector<std::string>(argv + 1, argv + argc);
    this->argv_i = argv_.begin();
    this->core();
}


DEploidIO::~DEploidIO() {
    if (this->isCopied()){
        return;
    }
    if ( this->excludedMarkers != NULL ) {
        delete this->excludedMarkers;
    }
    if ( this->vcfReaderPtr_ != NULL ) {
        delete this->vcfReaderPtr_;
    }
    if ( this->panel != NULL ) {
        delete panel;
    }
}


void DEploidIO::core() {
    if ( argv_.size() == 0 ) {
        this->set_help(true);
        return;
    }

    this->reInit(); // Reset to defalt values before parsing
    this->parse();

    if ( this->help() || version() ) {
        return;
    }

    this->checkInput();
    this->finalize();
}


void DEploidIO::init() {
    this->setDoPrintLassoPanel(false);
    this->setIsCopied(false);
    this->setDoExportRecombProb(false);

    this->setCompressVcf(false);
    this->setInitialPropWasGiven(false);
    this->setInitialHapWasGiven(false);
    this->initialProp.clear();
    this->setPleaseCheckInitialP(true);
    this->setExcludeSites( false );
    this->excludedMarkers = NULL;
    this->panel = NULL;
    this->set_help(false);
    this->setVersion(false);
    this->setUsePanel(true);
    this->prefix_ = "pf3k-dEploid";
    this->setKStrainWasSetByHap(false);
    this->setKStrainWasSetByProp(false);

    this->setDoUpdateProp( true );
    this->setDoUpdatePair( true );
    this->setDoUpdateSingle( true );
    this->setDoExportPostProb( false );
    this->setDoLsPainting( false );
    this->setDoIbdPainting( false );
    this->setDoIbdViterbiPainting(false);
    this->setUseIBD( false );
    this->setUseIbdOnly(false);
    this->setUseLasso( false );
    this->setUseBestPractice(false);
    this->setInferBestPracticeP(true);
    this->setInferBestPracticeHap(true);
    this->setDoExportSwitchMissCopy ( true );
    this->setDoAllowInbreeding( false );
    this->useConstRecomb_ = false;
    this->setForbidCopyFromSame( false );
    this->constRecombProb_ = 1.0;
    this->averageCentimorganDistance_ = 15000.0;
    this->setScalingFactor(100.0);
    this->setParameterG(20.0);
    this->setIBDSigma(20.0);
    this->setUseVcf(false);
    this->setUseVcfSample(false);
    this->setExtractPlafFromVcf(false);
    this->vcfSampleName_ = "";
    this->vcfReaderPtr_ = NULL;
    this->setDoExportVcf(false);
    this->setDoComputeLLK(false);
    this->setVqslod(8.0);
    this->setLassoMaxNumPanel(100);

    this->kStrain_.init(5);  // From DEploid-Lasso, set default K to 4.
    this->mcmcBurn_.init(0.5);
    this->mcmcMachineryRate_.init(5);
    this->missCopyProb_.init(0.01);
    this->nMcmcSample_.init(800);
    this->precision_.init(8);
    this->randomSeed_.init((unsigned)0);
    this->parameterSigma_.init(5.0);


    #ifdef COMPILEDATE
        compileTime_ = COMPILEDATE;
    #else
        compileTime_ = "";
    #endif

    #ifdef DEPLOIDVERSION
        dEploidGitVersion_ = DEPLOIDVERSION;
    #else
        dEploidGitVersion_ = "";
    #endif

    #ifdef DEPLOIDVERSION
        lassoGitVersion_ = LASSOVERSION;
    #else
        lassoGitVersion_ = "";
    #endif


}


void DEploidIO::setBestPracticeParameters() {
    this->kStrain_.setBest(4);
    this->nMcmcSample_.setBest(500);
    this->randomSeed_.setBest(1);
    this->parameterSigma_.setBest(1.6);
    this->mcmcBurn_.setBest(.67);
    this->mcmcMachineryRate_.setBest(8);
}

void DEploidIO::getTime( bool isStartingTime ) {
    time_t now = time(0);
    // convert now to string form
    char* dt = ctime(&now);
    if ( isStartingTime ) {
        startingTime_ = string(dt);
    } else {
        endTime_ = string(dt);
    }
}


void DEploidIO::reInit() {
    this->init();
    this->refFileName_.clear();
    this->altFileName_.clear();
    this->plafFileName_.clear();
    this->panelFileName_.clear();
    this->excludeFileName_.clear();
    this->getTime(true);
}


void DEploidIO::finalize() {
    if (useBestPractice()){
        this->setBestPracticeParameters();
    }

    if ( this->doIbdPainting() || this->doComputeLLK() || this->doIbdViterbiPainting() ) {
        if (!initialPropWasGiven()) {
            throw InitialPropUngiven("");
        }
    }

    if ( this->useIBD() && this->kStrain_.getValue() == 1) {
        throw InvalidK();
    }

    if ( this->compressVcf() && !this->doExportVcf() ) {
        throw VcfOutUnSpecified("");
    }

    if ( !this->randomSeed_.useDefault() ) {
        this->randomSeed_.init((unsigned)(time(0)));
    }

    if ( this->excludeSites() ) {
        excludedMarkers = new ExcludeMarker();
        excludedMarkers->readFromFile(excludeFileName_.c_str());
    }

    if ( useVcf() ) { // read vcf files, and parse it to refCount and altCount
        this->vcfReaderPtr_ = new VcfReader (vcfFileName_, vcfSampleName_,
            extractPlafFromVcf_);
        if ( this->excludeSites() ) {
            this->vcfReaderPtr_->findAndKeepMarkers (excludedMarkers);
        }
        this->vcfReaderPtr_->finalize(); // Finalize after remove variantlines
        this->refCount_ = this->vcfReaderPtr_->refCount;
        this->altCount_ = this->vcfReaderPtr_->altCount;
    } else {
        TxtReader ref;
        ref.readFromFile(refFileName_.c_str());
        if ( this->excludeSites() ) {
            ref.findAndKeepMarkers( excludedMarkers );
        }
        this->refCount_ = ref.info_;

        TxtReader alt;
        alt.readFromFile(altFileName_.c_str());
        if ( this->excludeSites() ) {
            alt.findAndKeepMarkers( excludedMarkers );
        }
        this->altCount_ = alt.info_;
    }

    this->nLoci_ = refCount_.size();

    if ( this->nLoci_ != this->altCount_.size() ) {
        throw LociNumberUnequal( this->altFileName_ );
    }

    if (extractPlafFromVcf()){
        this->plaf_ = this->vcfReaderPtr_->plaf;
        this->chrom_ = this->vcfReaderPtr_->chrom_;
        this->position_ = this->vcfReaderPtr_->position_;
        this->indexOfChromStarts_ = this->vcfReaderPtr_->indexOfChromStarts_;
    } else {
        TxtReader plaf;
        plaf.readFromFile(plafFileName_.c_str());
        if ( this->excludeSites() ) {
            plaf.findAndKeepMarkers( excludedMarkers );
        }
        this->plaf_ = plaf.info_;
        this->chrom_ = plaf.chrom_;
        this->position_ = plaf.position_;
        this->indexOfChromStarts_ = plaf.indexOfChromStarts_;
        if ( this->nLoci_ != this->plaf_.size() ) {
            throw LociNumberUnequal( this->plafFileName_ );
        }
    }

    (void)removeFilesWithSameName();

    this->readPanel();
    IBDpathChangeAt = vector <double>(this->nLoci());
    finalIBDpathChangeAt = vector <double>(this->nLoci());

    siteOfTwoSwitchOne = vector <double>(this->nLoci());
    siteOfTwoMissCopyOne = vector <double>(this->nLoci());
    siteOfTwoSwitchTwo = vector <double>(this->nLoci());
    siteOfTwoMissCopyTwo = vector <double>(this->nLoci());
    siteOfOneSwitchOne = vector <double>(this->nLoci());
    siteOfOneMissCopyOne = vector <double>(this->nLoci());

    finalSiteOfTwoSwitchOne = vector <double>(this->nLoci());
    finalSiteOfTwoMissCopyOne = vector <double>(this->nLoci());
    finalSiteOfTwoSwitchTwo = vector <double>(this->nLoci());
    finalSiteOfTwoMissCopyTwo = vector <double>(this->nLoci());
    finalSiteOfOneSwitchOne = vector <double>(this->nLoci());
    finalSiteOfOneMissCopyOne = vector <double>(this->nLoci());
}


void DEploidIO::removeFilesWithSameName() {
    //strExportProp = this->prefix_ + ".prop";
    //strExportLLK = this->prefix_ + ".llk";
    //strExportHap = this->prefix_ + ".hap";

    //if (this->useIbdOnly()) {
        //strIbdExportProp = this->prefix_ + ".prop";
        //strIbdExportLLK = this->prefix_ + ".llk";
        //strIbdExportHap = this->prefix_ + ".hap";
    //} else {
        //strIbdExportProp = this->prefix_ + ".ibd.prop";
        //strIbdExportLLK = this->prefix_ + ".ibd.llk";
        //strIbdExportHap = this->prefix_ + ".ibd.hap";
        strIbdExportProbs = this->prefix_ + ".ibd.probs";
        strIbdExportViterbi = this->prefix_ + ".viterbi";
    //}


    strExportLog =  this->prefix_ + ((this->doLsPainting()) ? ".painting":"") + ".log";
    strExportRecombProb = this->prefix_ + ".recomb";

    strExportExtra = this->prefix_ + ".extra";

    if ( this->doLsPainting() == false ) {
        if (this->useIBD()) {
            //remove(strIbdExportProp.c_str());
            //remove(strIbdExportLLK.c_str());
            //remove(strIbdExportHap.c_str());
        }
        //remove(strExportLLK.c_str());
        //remove(strExportHap.c_str());
        //remove(strExportVcf.c_str());
        //remove(strExportProp.c_str());
        remove(strExportExtra.c_str());
        remove(strIbdExportProbs.c_str());
        remove(strIbdExportViterbi.c_str());
    }

    if (this->doLsPainting() || this->doExportPostProb() ) {
        if (this->useIBD()) {
            strIbdExportSingleFwdProbPrefix = this->prefix_ + ".ibd.single";
            for ( size_t i = 0; i < this->kStrain_.getValue() ; i++ ) {
                string tmpStrExportSingleFwdProb = strIbdExportSingleFwdProbPrefix + to_string(i);
                remove(tmpStrExportSingleFwdProb.c_str());
            }
            strIbdExportPairFwdProb = this->prefix_ + ".ibd.pair";
            remove(strIbdExportPairFwdProb.c_str());
        }
        strExportSingleFwdProbPrefix = this->prefix_ + ".single";
        for ( size_t i = 0; i < this->kStrain_.getValue() ; i++ ) {
            string tmpStrExportSingleFwdProb = strExportSingleFwdProbPrefix + to_string(i);
            remove(tmpStrExportSingleFwdProb.c_str());
        }
        strExportPairFwdProb = this->prefix_ + ".pair";
        remove(strExportPairFwdProb.c_str());
    }
    remove(strExportLog.c_str());
    remove(strExportRecombProb.c_str());

}


void DEploidIO::parse () {
    do {
        if (*argv_i == "-ref") {
            if ( this->useVcf() ) {
                throw ( FlagsConflict((*argv_i) , "-vcf") );
            }
            this->readNextStringto ( this->refFileName_ ) ;
        } else if (*argv_i == "-alt") {
            if ( this->useVcf() ) {
                throw ( FlagsConflict((*argv_i) , "-vcf") );
            }
            this->readNextStringto ( this->altFileName_ ) ;
        } else if (*argv_i == "-vcf") {
            if ( this->refFileName_.size() > 0 || this->altFileName_.size() > 0 ) {
                throw ( FlagsConflict((*argv_i) , "-ref or -alt") );
            }
            this->setUseVcf(true);
            this->readNextStringto ( this->vcfFileName_ ) ;
        } else if (*argv_i == "-sample") {
            if (useVcf() == false){
                throw(FlagsOrderIncorrect((*argv_i) , "-vcf") );
            }
            this->setUseVcfSample(true);
            this->readNextStringto ( this->vcfSampleName_ ) ;
        } else if (*argv_i == "-plafFromVcf") {
            if (useVcf() == false){
                throw(FlagsOrderIncorrect((*argv_i) , "-vcf") );
            }
            this->setExtractPlafFromVcf(true);
        } else if (*argv_i == "-vcfOut") {
            this->setDoExportVcf (true);
        } else if (*argv_i == "-plaf") {
            if ( this->extractPlafFromVcf() ) {
                throw ( FlagsConflict((*argv_i) , "-plafFromVcf") );
            }
            this->readNextStringto ( this->plafFileName_ ) ;
        } else if (*argv_i == "-panel") {
            if ( this->usePanel() == false ) {
                throw ( FlagsConflict((*argv_i) , "-noPanel") );
            }
            this->readNextStringto ( this->panelFileName_ ) ;
        } else if (*argv_i == "-noPanel") {
            if ( usePanel() && this->panelFileName_.size() > 0 ) {
                throw ( FlagsConflict((*argv_i) , "-panel") );
            }
            if ( doExportPostProb() ) {
                throw ( FlagsConflict((*argv_i) , "-exportPostProb") );
            }
            if ( doAllowInbreeding() ) {
                throw ( FlagsConflict((*argv_i) , "-inbreeding") );
            }
            if ( useBestPractice() ) {
                throw ( FlagsConflict((*argv_i) , "-best") );
            }
            this->setUsePanel(false);
            this->setDoExportSwitchMissCopy ( false );
        } else if (*argv_i == "-exclude") {
            this->setExcludeSites( true );
            this->readNextStringto ( this->excludeFileName_ ) ;
        } else if (*argv_i == "-o") {
            this->readNextStringto ( this->prefix_ ) ;
        } else if ( *argv_i == "-p" ) {
            //this->precision_ = readNextInput<size_t>() ;
            this->precision_.setUserDefined(readNextInput<size_t>());
        } else if ( *argv_i == "-k" ) {
            this->kStrain_.setUserDefined(readNextInput<size_t>());

            if ( this->kStrainWasSetByHap() && this->kStrain_.getValue() != this->initialHap[0].size() ) {
                string hint = string(" k = ") + to_string(this->kStrain_.getValue()) + ", " + this->initialHapFileName_ + " suggests otherwise";
                throw NumOfPropNotMatchNumStrain(hint);
            }

            if ( this->initialPropWasGiven() && this->kStrain_.getValue() != initialProp.size() ) {
                string hint = string(" k = ") + to_string(kStrain_.getValue()) + ", flag -initialP suggests otherwise";;
                throw NumOfPropNotMatchNumStrain(hint);
            }

        } else if ( *argv_i == "-nSample" ) {
            this->nMcmcSample_.setUserDefined(readNextInput<size_t>());
        } else if ( *argv_i == "-burn" ) {
            this->mcmcBurn_.setUserDefined(readNextInput<double>());
            if ( this->mcmcBurn_.getValue() < 0 || this->mcmcBurn_.getValue() > 1) {
                throw ( OutOfRange ("-burn", *argv_i) );
            }
        } else if ( *argv_i == "-miss" ) {
            this->missCopyProb_.setUserDefined(readNextInput<double>());
            if ( this->missCopyProb_.getValue() < 0 || this->missCopyProb_.getValue() > 1) {
                throw ( OutOfRange ("-miss", *argv_i) );
            }
        } else if ( *argv_i == "-c" ) {
            this->scalingFactor_ = readNextInput<double>() ;
        } else if ( *argv_i == "-G" ) {
            this->setParameterG(readNextInput<double>());
        } else if ( *argv_i == "-vqslod" ) {
            this->setVqslod(readNextInput<double>());
        } else if ( *argv_i == "-sigma" ) {
            this->parameterSigma_.setUserDefined(readNextInput<double>());
        } else if ( *argv_i == "-ibdSigma" ) {
            this->setIBDSigma(readNextInput<double>());
        } else if ( *argv_i == "-lassoMaxPanel" ) {
            this->setLassoMaxNumPanel(readNextInput<size_t>());
        } else if ( *argv_i == "-recomb" ) {
            this->constRecombProb_ = readNextInput<double>();
            this->useConstRecomb_ = true;
            if ( this->constRecombProb_ < 0 || this->constRecombProb_ > 1) {
                throw ( OutOfRange ("-recomb", *argv_i) );
            }
        } else if ( *argv_i == "-printRecomb" ) {
            this->setDoExportRecombProb( true );
        } else if ( *argv_i == "-forbidSame" ) {
            this->setForbidCopyFromSame( true );
        } else if ( *argv_i == "-rate" ) {
            this->mcmcMachineryRate_.setUserDefined(readNextInput<size_t>());
        } else if ( *argv_i == "-forbidUpdateProp" ) {
            this->setDoUpdateProp( false );
        } else if ( *argv_i == "-forbidUpdateSingle" ) {
            this->setDoUpdateSingle( false );
        } else if ( *argv_i == "-forbidUpdatePair" ) {
            this->setDoUpdatePair( false );
        } else if ( *argv_i == "-inbreeding" ) {
            if ( this->usePanel() == false ) {
                throw ( FlagsConflict((*argv_i) , "-noPanel") );
            }
            this->setDoAllowInbreeding( true );
            this->setDoExportPostProb( true );
        } else if ( *argv_i == "-exportPostProb" ) {
            if ( useBestPractice() ) {
                throw ( FlagsConflict((*argv_i) , "-best") );
            }
            if ( this->usePanel() == false ) {
                throw ( FlagsConflict((*argv_i) , "-noPanel") );
            }
            this->setDoExportPostProb( true );
        } else if ( *argv_i == "-painting" ) {
            if ( this->usePanel() == false ) {
                throw ( FlagsConflict((*argv_i) , "-noPanel") );
            }
            if ( this->initialHapWasGiven() == true ) {
                throw ( FlagsConflict((*argv_i) , "-initialHap") );
            }
            this->readNextStringto ( this->initialHapFileName_ ) ;
            this->setDoLsPainting( true );
            this->readInitialHaps();
        } else if ( *argv_i == "-ibd" ) {
            if ( useBestPractice() ) {
                throw ( FlagsConflict((*argv_i) , "-best") );
            }
            if ( useLasso() == true ) {
                throw ( FlagsConflict((*argv_i) , "-lasso") );
            }
            this->setUseIBD(true);
        } else if (*argv_i == "-best") {
            if ( doExportPostProb() ) {
                throw ( FlagsConflict((*argv_i) , "-exportPostProb") );
            }
            if ( usePanel() == false ) {
                throw ( FlagsConflict((*argv_i) , "-noPanel") );
            }
            if ( useIBD() == true ) {
                throw ( FlagsConflict((*argv_i) , "-ibd") );
            }
            if ( useLasso() == true ) {
                throw ( FlagsConflict((*argv_i) , "-lasso") );
            }
            this->setUseBestPractice(true);
        } else if (*argv_i == "-bestKonly") {
            this->setInferBestPracticeP(false);
            this->setInferBestPracticeHap(false);
        } else if (*argv_i == "-bestPonly") {
            this->setInferBestPracticeHap(false);
        } else if ( *argv_i == "-ibdonly" ) {
            if ( useBestPractice() ) {
                throw ( FlagsConflict((*argv_i) , "-best") );
            }
            this->setUseIBD(true);
            this->setUseIbdOnly(true);
        } else if ( *argv_i == "-lasso" ) {
            if ( useIBD() == true ) {
                throw ( FlagsConflict((*argv_i) , "-ibd") );
            }
            if ( useBestPractice() ) {
                throw ( FlagsConflict((*argv_i) , "-best") );
            }
            if ( doExportPostProb() ) {
                throw ( FlagsConflict((*argv_i) , "-exportPostProb") );
            }
            this->setUseLasso(true);
            this->setDoUpdateProp(false);
        } else if ( *argv_i == "-writePanel" ) {
            this->setDoPrintLassoPanel(true);
        } else if ( *argv_i == "-computeLLK" ) {
            this->setDoComputeLLK( true );
        } else if ( *argv_i == "-ibdPainting" ) {
            if ( useBestPractice() ) {
                throw ( FlagsConflict((*argv_i) , "-best") );
            }
            this->setDoIbdPainting( true );
        } else if ( *argv_i == "-ibdViterbi" ) {
            this->setDoIbdViterbiPainting( true );
        } else if ( *argv_i == "-skipCheckingInitialP" ) {
            this->setPleaseCheckInitialP(false);
        } else if ( *argv_i == "-initialP" ) {
            this->readInitialProportions();
            this->setInitialPropWasGiven( true );

            // If the k was set manually, check
            if ( this->kStrain_.useUserDefined() && this->kStrain_.getValue() != initialProp.size() ) {
                string hint = string(" k = ") + to_string(kStrain_.getValue());
                throw NumOfPropNotMatchNumStrain(hint);
            }

            // If the k was set by initial Hap, check
            if ( this->kStrainWasSetByHap() && this->kStrain_.getValue() != this->initialProp.size() ) {
                string hint = string(" k = ") + to_string(this->kStrain_.getValue()) + ", " + this->initialHapFileName_ + " suggests otherwise";
                throw NumOfPropNotMatchNumStrain(hint);
            }

            this->kStrain_.init(this->initialProp.size());
        } else if ( *argv_i == "-initialHap" ) {
            if ( this->doLsPainting() == true ) {
                throw ( FlagsConflict((*argv_i) , "-painting") );
            }
            this->readNextStringto ( this->initialHapFileName_ ) ;
            this->setInitialHapWasGiven(true);
            this->readInitialHaps();
        } else if ( *argv_i == "-seed") {
            this->randomSeed_.setUserDefined(readNextInput<size_t>());
        } else if ( *argv_i == "-z" ) {
            this->setCompressVcf(true);
        } else if ( *argv_i == "-h" || *argv_i == "-help") {
            this->set_help(true);
        } else if ( *argv_i == "-v" || *argv_i == "-version") {
            this->setVersion(true);
        } else {
            throw ( UnknowArg((*argv_i)) );
        }
    } while ( ++argv_i != argv_.end());
}


void DEploidIO::checkInput() {

    if ( this->refFileName_.size() == 0 && this->useVcf() == false ) {
        throw FileNameMissing ( "Ref count" );}
    if ( this->altFileName_.size() == 0 && this->useVcf() == false ) {
        throw FileNameMissing ( "Alt count" );}
    if ( this->plafFileName_.size() == 0 && extractPlafFromVcf() == false ) {
        throw FileNameMissing ( "PLAF" );}
    if ( usePanel() && this->panelFileName_.size() == 0 && !this->doIbdPainting() && !this->doComputeLLK() ) {
        throw FileNameMissing ( "Reference panel" );}
    if ( this->initialPropWasGiven() && ( abs(sumOfVec(initialProp) - 1.0) > 0.00001 ) && this->pleaseCheckInitialP() ) {
        throw SumOfPropNotOne ( to_string(sumOfVec(initialProp)) );}
    if ( this->initialPropWasGiven() ) {
        if ( this->kStrain_.useUserDefined() == true ) {
        } else {
            // set k strain by proportion length
        }
    }
    if (this->useBestPractice() && (!this->usePanel())){
        throw FlagsConflict("-best" , string("-noPanel. Reference panel is") +
            string("required for using best-practices."));
    }
}


void DEploidIO::readInitialProportions() {
    string tmpFlag = *argv_i;
    ++argv_i;
    if (argv_i == argv_.end() || (*argv_i)[0] == '-' ) {
            throw NotEnoughArg (tmpFlag);
    }

    do {
        try {
            double tmp = convert<double>(*argv_i);
            this->initialProp.push_back(tmp);
        } catch (const WrongType &e) {
            --argv_i;
            break;
        }
        ++argv_i;
    } while ( argv_i != argv_.end() && (*argv_i)[0] != '-' );
    --argv_i;

    return;
}


void DEploidIO::readNextStringto( string &readto ) {
    string tmpFlag = *argv_i;
    ++argv_i;
    if (argv_i == argv_.end() || (*argv_i)[0] == '-' ) {
        throw NotEnoughArg(tmpFlag);
    }
    readto = *argv_i;
}



std::ostream& operator<< (std::ostream& stream, const DEploidIO& dEploidIO) {
  for (std::string arg : dEploidIO.argv_) stream << " " << arg;
  return stream;
}


void DEploidIO::readInitialHaps() {
    assert( this->initialHap.size() == 0 );
    InitialHaplotypes initialHapToBeRead;
    initialHapToBeRead.readFromFile(this->initialHapFileName_.c_str());

    assert (this->initialHap.size() == 0 );
    this->initialHap = initialHapToBeRead.content_;

    if ( this->kStrain_.useUserDefined() && this->kStrain_.getValue()!= initialHapToBeRead.truePanelSize() ) {
        string hint = string(" k = ") + to_string(this->kStrain_.getValue()) + ", " + this->initialHapFileName_ + " suggests otherwise";
        throw NumOfPropNotMatchNumStrain(hint);
    }

    if ( this->kStrainWasSetByProp() && this->kStrain_.getValue() != initialHapToBeRead.truePanelSize() ) {
        string hint = string(" k = ") + to_string(kStrain_.getValue()) + ", flag -initialP suggests otherwise";;
        throw NumOfPropNotMatchNumStrain(hint);
    }

    this->kStrain_.init(initialHapToBeRead.truePanelSize());
    this->setKStrainWasSetByHap(true);
}


vector <double> DEploidIO::computeExpectedWsafFromInitialHap() {
    // Make this a separate function
    // calculate expected wsaf
    vector <double> expectedWsaf (this->initialHap.size(), 0.0);
    for ( size_t i = 0; i < this->initialHap.size(); i++ ) {
        assert( kStrain_.getValue() == this->initialHap[i].size() );
        for ( size_t k = 0; k < this->kStrain_.getValue(); k++) {
            expectedWsaf[i] += this->initialHap[i][k] * finalProp[k];
        }
        assert ( expectedWsaf[i] >= 0 );
        //assert ( expectedWsaf[i] <= 1.0 );
    }
    return expectedWsaf;
}


void DEploidIO::computeLLKfromInitialHap() {
    for ( auto const& value: this->initialProp ) {
        this->finalProp.push_back(value);
    }

    vector <double> expectedWsaf = computeExpectedWsafFromInitialHap();
    if (expectedWsaf.size() != this->refCount_.size()) {
        throw LociNumberUnequal("Hap length differs from data!");
    }
    vector <double> llk = calcLLKs ( this->refCount_, this->altCount_, expectedWsaf, 0, expectedWsaf.size(), this->scalingFactor());
    this->llkFromInitialHap_ = sumOfVec(llk);
}




void DEploidIO::readPanel() {
    if ( this->usePanel() == false ) {
        return;
    }
    if ( this->doIbdPainting() || this->doComputeLLK() ) {
        return;
    }

    panel = new Panel();
    panel->readFromFile(this->panelFileName_.c_str());
    if ( this->excludeSites() ) {
        panel->findAndKeepMarkers( this->excludedMarkers );
    }

    panel->computeRecombProbs( this->averageCentimorganDistance(), this->parameterG(), this->useConstRecomb(), this->constRecombProb(), this->forbidCopyFromSame() );
    panel->checkForExceptions( this->nLoci(), this->panelFileName_ );
}


DEploidIO::DEploidIO(const DEploidIO &cpFrom) {
    this->setIsCopied(true);
    this->setDoExportRecombProb(cpFrom.doExportRecombProb());
    this->setCompressVcf(cpFrom.compressVcf());
    this->setInitialPropWasGiven(cpFrom.initialPropWasGiven());
    this->setInitialHapWasGiven(cpFrom.initialHapWasGiven());
    this->initialProp = vector <double> (cpFrom.initialProp.begin(),
                                         cpFrom.initialProp.end());
    this->setPleaseCheckInitialP(cpFrom.pleaseCheckInitialP());
    this->setExcludeSites(cpFrom.excludeSites());
    this->excludedMarkers = cpFrom.excludedMarkers;
    this->panel = cpFrom.panel;
    this->set_help(cpFrom.help());
    this->setVersion(cpFrom.version());
    this->setUsePanel(cpFrom.usePanel());
    this->prefix_ = cpFrom.prefix_;
    this->setKStrainWasSetByHap(cpFrom.kStrainWasSetByHap());
    this->setKStrainWasSetByProp(cpFrom.kStrainWasSetByProp());

    this->setDoUpdateProp(cpFrom.doUpdateProp());
    this->setDoUpdatePair(cpFrom.doUpdatePair());
    this->setDoUpdateSingle(cpFrom.doUpdateSingle());
    this->setDoExportPostProb(cpFrom.doExportPostProb());
    this->setDoLsPainting(cpFrom.doLsPainting());
    this->setDoIbdPainting(cpFrom.doIbdPainting());
    this->setUseIBD(cpFrom.useIBD());
    this->setUseIbdOnly(cpFrom.useIbdOnly());
    this->setUseLasso(cpFrom.useLasso());
    this->setDoExportSwitchMissCopy(cpFrom.doExportSwitchMissCopy());
    this->setDoAllowInbreeding(cpFrom.doAllowInbreeding());
    this->useConstRecomb_ = cpFrom.useConstRecomb();
    this->setForbidCopyFromSame(cpFrom.forbidCopyFromSame());
    this->constRecombProb_ = cpFrom.constRecombProb();
    this->averageCentimorganDistance_ = cpFrom.averageCentimorganDistance();
    this->setScalingFactor(cpFrom.scalingFactor());
    this->setParameterG(cpFrom.parameterG());
    this->setIBDSigma(cpFrom.ibdSigma());
    this->setUseVcf(cpFrom.useVcf());
    this->vcfReaderPtr_ = cpFrom.vcfReaderPtr_;
    this->setDoExportVcf(cpFrom.doExportVcf());
    this->setDoComputeLLK(cpFrom.doComputeLLK());
    this->setNLoci(cpFrom.nLoci());
    this->refCount_ = vector <double> (cpFrom.refCount_.begin(),
                                       cpFrom.refCount_.end());
    this->altCount_ = vector <double> (cpFrom.altCount_.begin(),
                                       cpFrom.altCount_.end());
    this->plaf_ = vector <double> (cpFrom.plaf_.begin(),
                                   cpFrom.plaf_.end());
    this->chrom_ = vector <string> (cpFrom.chrom_.begin(),
                                   cpFrom.chrom_.end());
    this->position_ = vector < vector <int> > (cpFrom.position_.begin(),
                                   cpFrom.position_.end());
    this->indexOfChromStarts_ = vector <size_t> (cpFrom.indexOfChromStarts_.begin(),
                                   cpFrom.indexOfChromStarts_.end());
    this->setVqslod(cpFrom.vqslod());
    this->setLassoMaxNumPanel(cpFrom.lassoMaxNumPanel());
    //this->strExportProp = cpFrom.strExportProp;
    //this->strExportLLK = cpFrom.strExportLLK;
    //this->strExportHap = cpFrom.strExportHap;
    //this->strIbdExportProp = cpFrom.strIbdExportProp;
    //this->strIbdExportLLK = cpFrom.strIbdExportLLK;
    //this->strIbdExportHap = cpFrom.strIbdExportHap;

    this->kStrain_.makeCopy(cpFrom.kStrain_);
    this->precision_.makeCopy(cpFrom.precision_);
    this->missCopyProb_ .makeCopy(cpFrom.missCopyProb_);
    this->randomSeed_.makeCopy(cpFrom.randomSeed_);
    this->nMcmcSample_.makeCopy(cpFrom.nMcmcSample_);
    this->parameterSigma_.makeCopy(cpFrom.parameterSigma_);
    this->mcmcBurn_.makeCopy(cpFrom.mcmcBurn_);
    this->mcmcMachineryRate_.makeCopy(cpFrom.mcmcMachineryRate_);
}


void DEploidIO::getIBDprobsIntegrated(vector < vector <double> > &prob) {
    if (prob.size() !=  this->nLoci()) {
        throw InvalidInput("Invalid probabilities! Check size!");
    }

    assert(this->ibdProbsIntegrated.size() == 0);

    for (size_t i = 0; i < prob[0].size(); i++) {
        this->ibdProbsIntegrated.push_back(0.0);
    }

    for ( size_t siteIndex = 0; siteIndex < this->nLoci(); siteIndex++ ) {
        for (size_t i = 0; i < prob[siteIndex].size(); i++) {
            this->ibdProbsIntegrated[i] += prob[siteIndex][i];
        }
    }
    normalizeBySum(this->ibdProbsIntegrated);
}


void DEploidIO::computeEffectiveKstrain(vector <double> proportion) {
    double tmpSumSq = 0.0;
    for (double p : proportion) {
        tmpSumSq += p * p;
    }
    this->effectiveKstrain_ = 1.0 / tmpSumSq;
}


void DEploidIO::computeInferredKstrain(vector <double> proportion) {
    this->inferredKstrain_ = 0;
    for (double p : proportion) {
        if ( p > 0.01 ) {
            this->inferredKstrain_ += 1;
        }
    }
}


void DEploidIO::computeAdjustedEffectiveKstrain() {
    this->adjustedEffectiveKstrain_ = this->effectiveKstrain_;
    if ( (this->inferredKstrain_ == 2) & (ibdProbsIntegrated.size() == 2)) {
        if ( this->ibdProbsIntegrated[1] > 0.95 ) {
            this->adjustedEffectiveKstrain_ = 1;
        }
    }
}


vector <double> DEploidIO::lassoComputeObsWsaf(size_t segmentStartIndex, size_t nLoci) {
    vector <double> ret(nLoci, 0.0);
    for ( size_t i = 0; i < nLoci; i++) {
        size_t idx = i + segmentStartIndex;
        ret[i] = this->altCount_[idx] / (this->refCount_[idx] + this->altCount_[idx] + 0.00000000000001);
    }
    return ret;
}


vector < vector <double> > DEploidIO::lassoSubsetPanel(size_t segmentStartIndex, size_t nLoci) {
    vector < vector <double> > ret;
    for ( size_t i = 0; i < nLoci; i++) {
        size_t idx = i + segmentStartIndex;
        ret.push_back(this->panel->content_[idx]);
    }
    return ret;
}


void DEploidIO::dEploidLasso() {
    for ( size_t chromi = 0 ; chromi < this->indexOfChromStarts_.size(); chromi++ ) {
        size_t start = this->indexOfChromStarts_[chromi];
        size_t length = this->position_[chromi].size();
        size_t end = start + length;
        vector <double> wsaf = lassoComputeObsWsaf(start, length);
        vector < vector <double> > tmpPanel = lassoSubsetPanel(start, length);
        DEploidLASSO dummy(tmpPanel, wsaf, 250);

        vector <string> newHeader;
        for (size_t i = 0; i < min(dummy.choiceIdx.size(), lassoMaxNumPanel()); i++) {
            newHeader.push_back(panel->header_[dummy.choiceIdx[i]]);
        }
        newHeader.push_back("3d7");

        vector < vector <double> > newPanel;
        for (size_t i = 0; i < dummy.reducedPanel.size(); i++) {
            vector <double> tmpRow;
            for (size_t j = 0; j < min(dummy.choiceIdx.size(), lassoMaxNumPanel()); j++) {
                tmpRow.push_back(dummy.reducedPanel[i][j]);
            }
            tmpRow.push_back(static_cast<int>(0));
            newPanel.push_back(tmpRow);
        }

        Panel * tmp = new Panel(vecFromTo(this->panel->pRec_, start, end),
                            vecFromTo(this->panel->pRecEachHap_, start, end),
                            vecFromTo(this->panel->pNoRec_, start, end),
                            vecFromTo(this->panel->pRecRec_, start, end),
                            vecFromTo(this->panel->pRecNoRec_, start, end),
                            vecFromTo(this->panel->pNoRecNoRec_, start, end),
                            newPanel,
                            this->panel->header_);
        lassoPanels.push_back(tmp);
        lassoPlafs.push_back(vecFromTo(plaf_, start, end));
        lassoRefCount.push_back(vecFromTo(refCount_, start, end));
        lassoAltCount.push_back(vecFromTo(altCount_, start, end));
        if (this->doPrintLassoPanel()) {
            this->writePanel(tmp, chromi, newHeader);
        }
    }
    this->finalProp.clear();
    for ( auto const& value: this->initialProp ) {
        this->finalProp.push_back(value);
    }
}


void DEploidIO::trimVec(vector <double> &vec, vector <size_t> &idx) {
    vector <double> ret;
    for (auto const& value : idx){
        ret.push_back(vec[value]);
    }
    //return ret;
    vec.clear();
    for (auto const& value : ret){
        vec.push_back(value);
    }
}


void DEploidIO::trimming(vector <size_t> & trimmingCriteria) {
    this->trimVec(this->refCount_, trimmingCriteria);
    this->trimVec(this->altCount_, trimmingCriteria);
    this->trimVec(this->plaf_, trimmingCriteria);

    this->setNLoci(this->plaf_.size());

    vector <string> oldChrom = vector <string> (chrom_.begin(), chrom_.end());
    this->chrom_.clear();

    vector < vector < int > > oldposition = this->position_;
    this->position_.clear();

    for (size_t chromI = 0; chromI < oldChrom.size(); chromI++) {
        size_t hapIndex = indexOfChromStarts_[chromI];
        vector <int> newTrimmedPos;
        for (size_t posI = 0; posI < oldposition[chromI].size(); posI++) {
            if (std::find(trimmingCriteria.begin(),trimmingCriteria.end(), hapIndex)
                    != trimmingCriteria.end()){
                if (newTrimmedPos.size() == 0) {
                    this->chrom_.push_back(oldChrom[chromI]);
                }
                newTrimmedPos.push_back(oldposition[chromI][posI]);
            }

            hapIndex++;
        }
        this->position_.push_back(newTrimmedPos);
    }

    this->indexOfChromStarts_.clear();
    assert(indexOfChromStarts_.size() == 0);
    this->indexOfChromStarts_.push_back((size_t)0);
    for (size_t tmpChrom = 0;
            indexOfChromStarts_.size() < this->chrom_.size(); tmpChrom++ ) {
        indexOfChromStarts_.push_back(
            indexOfChromStarts_.back()+this->position_[tmpChrom].size());
    }
    assert(indexOfChromStarts_.size() == this->chrom_.size());
}



void DEploidIO::trimmingHalf(vector <size_t> & trimmingCriteria) {
    this->trimVec(this->refCount_, trimmingCriteria);
    this->trimVec(this->altCount_, trimmingCriteria);
    this->trimVec(this->plaf_, trimmingCriteria);

    this->setNLoci(this->plaf_.size());

    vector <string> oldChrom = vector <string> (chrom_.begin(), chrom_.end());
    this->chrom_.clear();

    vector < vector < int > > oldposition = this->position_;
    this->position_.clear();

    for (size_t chromI = 0; chromI < oldChrom.size(); chromI++) {
        if (chromI > 10) {
        //if (chromI%2 == 0) {
            size_t hapIndex = indexOfChromStarts_[chromI];
            vector <int> newTrimmedPos;
            for (size_t posI = 0; posI < oldposition[chromI].size(); posI++) {
                if (std::find(trimmingCriteria.begin(),trimmingCriteria.end(), hapIndex)
                        != trimmingCriteria.end()){
                    if (newTrimmedPos.size() == 0) {
                        this->chrom_.push_back(oldChrom[chromI]);
                    }
                    newTrimmedPos.push_back(oldposition[chromI][posI]);
                }

                hapIndex++;
            }
            this->position_.push_back(newTrimmedPos);
        }
    }

    this->indexOfChromStarts_.clear();
    assert(indexOfChromStarts_.size() == 0);
    this->indexOfChromStarts_.push_back((size_t)0);
    for (size_t tmpChrom = 0;
            indexOfChromStarts_.size() < this->chrom_.size(); tmpChrom++ ) {
        indexOfChromStarts_.push_back(
            indexOfChromStarts_.back()+this->position_[tmpChrom].size());
    }
    assert(indexOfChromStarts_.size() == this->chrom_.size());
}


void DEploidIO::computeObsWsaf() {
    assert(this->obsWsaf_.size() == 0);
    for ( size_t i = 0; i < this->nLoci(); i++) {
        this->obsWsaf_.push_back(this->altCount_[i] /
            (this->refCount_[i] + this->altCount_[i] + 0.00000000000001));
    }
    assert(this->obsWsaf_.size() == this->nLoci());
}


void DEploidIO::dEploidLassoTrimfirst() {  // This is trimming using VQSLOD
    if (vcfReaderPtr_ == NULL){
        panel = new Panel (*panel);
        panel->computeRecombProbs(this->averageCentimorganDistance(), this->parameterG(), true, 0.0000001, this->forbidCopyFromSame());
        this->dEploidLasso();
        this->setIsCopied(false);
        this->excludedMarkers = NULL;
        return;
    }

    this->vcfReaderPtr_->findLegitSnpsGivenVQSLOD(this->vqslod());
    this->trimming(this->vcfReaderPtr_->legitVqslodAt);

    panel = new Panel (*panel);
    panel->computeRecombProbs(this->averageCentimorganDistance(), this->parameterG(), true, 0.0000001, this->forbidCopyFromSame());
    panel->findAndKeepMarkersGivenIndex(this->vcfReaderPtr_->legitVqslodAt);
    this->dEploidLasso();
    this->setIsCopied(false);
    this->excludedMarkers = NULL;
    this->vcfReaderPtr_ = NULL;
}



void DEploidIO::dEploidLassoFullPanel() {
    // Filter SNPs first, restrict to wsaf > 0 don't work ...
    this->vcfReaderPtr_->findLegitSnpsGivenVQSLODHalf(this->vqslod());
    //this->vcfReaderPtr_->findLegitSnpsGivenVQSLOD(this->vqslod());
    //this->vcfReaderPtr_->findLegitSnpsGivenVQSLODandWsfGt0(this->vqslod());

    //this->trimming(this->vcfReaderPtr_->legitVqslodAt);
    this->trimmingHalf(this->vcfReaderPtr_->legitVqslodAt);
    this->computeObsWsaf();

    Panel tmpPanel(*panel);
    tmpPanel.computeRecombProbs(this->averageCentimorganDistance(), this->parameterG(), true, 0.0000001, this->forbidCopyFromSame());
    //tmpPanel.findAndKeepMarkersGivenIndex(this->vcfReaderPtr_->legitVqslodAt);
    tmpPanel.findAndKeepMarkersGivenIndexHalf(this->vcfReaderPtr_->legitVqslodAt);
    DEploidLASSO dummy(tmpPanel.content_, this->obsWsaf_, 250);

    for (size_t i = 0; i < dummy.choiceIdx.size(); i++) {
        dout << i << " " << dummy.devRatio[i]<<endl;
    }
    //size_t maxNumPanel = 10;
    // Use the first 10 strains in the panel
    vector <string> newHeader;
    for (size_t i = 0; i < min(dummy.choiceIdx.size(), lassoMaxNumPanel()); i++) {
        newHeader.push_back(panel->header_[dummy.choiceIdx[i]]);
    }
newHeader.push_back("3d7");
//cout <<newHeader.size()<<endl;

    vector < vector <double> > newPanel;
    for (size_t i = 0; i < dummy.reducedPanel.size(); i++) {
        vector <double> tmpRow;
        for (size_t j = 0; j < min(dummy.choiceIdx.size(), lassoMaxNumPanel()); j++) {
            tmpRow.push_back(dummy.reducedPanel[i][j]);
        }
tmpRow.push_back(static_cast<int>(0));
        newPanel.push_back(tmpRow);
    }

    panel = new Panel(tmpPanel.pRec_,
                      tmpPanel.pRecEachHap_,
                      tmpPanel.pNoRec_,
                      tmpPanel.pRecRec_,
                      tmpPanel.pRecNoRec_,
                      tmpPanel.pNoRecNoRec_,
                      newPanel,
                      newHeader);

    this->setIsCopied(false);
    this->excludedMarkers = NULL;
    this->vcfReaderPtr_ = NULL;
}


void DEploidIO::ibdTrimming() {
    if (vcfReaderPtr_ == NULL){
        panel = new Panel(*panel);
        this->setIsCopied(false);
        this->excludedMarkers = NULL;
        return;
    }
    // Filter SNPs first
    //cout << "here" <<endl;
    this->vcfReaderPtr_->findLegitSnpsGivenVQSLOD(this->vqslod());
    //cout << "stop here" <<endl;
    this->trimming(this->vcfReaderPtr_->legitVqslodAt);
    panel = new Panel(*panel);
    this->panel->findAndKeepMarkersGivenIndex(
                                        this->vcfReaderPtr_->legitVqslodAt);

    this->setIsCopied(false);
    this->excludedMarkers = NULL;
    this->vcfReaderPtr_ = NULL;
}
