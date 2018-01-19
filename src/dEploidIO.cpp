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
#include <iomanip>        // std::setw
#include "utility.hpp"    // normailize by sum
#include "updateHap.hpp"  // chromPainting
#include "dEploidIO.hpp"
#include "ibd.hpp"

DEploidIO::DEploidIO(){
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


DEploidIO::~DEploidIO(){
    if ( this->excludedMarkers != NULL ){
        delete this->excludedMarkers;
    }

    if ( this->vcfReaderPtr_ != NULL ){
        delete this->vcfReaderPtr_;
    }

    if ( this->panel != NULL ){
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

    if ( this->help() || version() ){
        return;
    }

    this->checkInput();
    this->finalize();
}


void DEploidIO::init() {
    this->setDoExportRecombProb(false);
    this->setrandomSeedWasGiven(false);
    this->setCompressVcf(false);
    this->setInitialPropWasGiven(false);
    this->setInitialHapWasGiven(false);
    this->initialProp.clear();
    this->setPleaseCheckInitialP(true);
    this->setExcludeSites( false );
    this->excludedMarkers = NULL;
    this->panel = NULL;
    this->set_seed( (unsigned)0 );
    this->set_help(false);
    this->setVersion(false);
    this->setUsePanel(true);
    this->precision_ = 8;
    this->prefix_ = "pf3k-dEploid";
    this->setKStrainWasManuallySet(false);
    this->setKStrainWasSetByHap(false);
    this->setKStrainWasSetByProp(false);
    this->setKstrain(5);
    this->nMcmcSample_ = 800;
    this->setDoUpdateProp( true );
    this->setDoUpdatePair( true );
    this->setDoUpdateSingle( true );
    this->setDoExportPostProb( false );
    this->setDoLsPainting( false );
    this->setDoIbdPainting( false );
    this->setUseIBD( false );
    this->setDoExportSwitchMissCopy ( true );
    this->setDoAllowInbreeding( false );
    this->mcmcBurn_ = 0.5;
    this->mcmcMachineryRate_ = 5;
    this->missCopyProb_ = 0.01;
    this->useConstRecomb_ = false;
    this->setForbidCopyFromSame( false );
    this->constRecombProb_ = 1.0;
    this->averageCentimorganDistance_ = 15000.0;
    this->setScalingFactor(100.0);
    this->setParameterG(20.0);
    this->setParameterSigma(5.0);
    this->setIBDSigma(20.0);
    this->setUseVcf(false);
    this->vcfReaderPtr_ = NULL;
    this->setDoExportVcf(false);
    this->setDoComputeLLK(false);

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
}


void DEploidIO::getTime( bool isStartingTime ){
    time_t now = time(0);
    // convert now to string form
    char* dt = ctime(&now);
    if ( isStartingTime ){
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


void DEploidIO::finalize(){
    if ( this->doIbdPainting() | this->doComputeLLK() ){
        if (!initialPropWasGiven()){
            throw InitialPropUngiven("");
        }
    }

    if ( this->useIBD() && this->kStrain() == 1){
        throw InvalidK();
    }

    if ( this->compressVcf() && !this->doExportVcf() ){
        throw VcfOutUnSpecified("");
    }

    if ( !this->randomSeedWasGiven_ ){
        this->set_seed( (unsigned)(time(0)) );
    }

    if ( this->excludeSites() ){
        excludedMarkers = new ExcludeMarker();
        excludedMarkers->readFromFile(excludeFileName_.c_str());
    }

    if ( useVcf() ){ // read vcf files, and parse it to refCount and altCount
        this->vcfReaderPtr_ = new VcfReader (vcfFileName_ );
        if ( this->excludeSites() ){
            this->vcfReaderPtr_->findAndKeepMarkers (excludedMarkers);
        }

        this->vcfReaderPtr_->finalize(); // Finalize after remove variantlines
        this->refCount_ = this->vcfReaderPtr_->refCount;
        this->altCount_ = this->vcfReaderPtr_->altCount;
    } else {
        TxtReader ref;
        ref.readFromFile(refFileName_.c_str());
        if ( this->excludeSites() ){
            ref.findAndKeepMarkers( excludedMarkers );
        }
        this->refCount_ = ref.info_;

        TxtReader alt;
        alt.readFromFile(altFileName_.c_str());
        if ( this->excludeSites() ){
            alt.findAndKeepMarkers( excludedMarkers );
        }
        this->altCount_ = alt.info_;
    }

    TxtReader plaf;
    plaf.readFromFile(plafFileName_.c_str());
    if ( this->excludeSites() ){
        plaf.findAndKeepMarkers( excludedMarkers );
    }
    this->plaf_ = plaf.info_;
    this->chrom_ = plaf.chrom_;
    this->position_ = plaf.position_;
    this->indexOfChromStarts_ = plaf.indexOfChromStarts_;

    this->nLoci_ = refCount_.size();

    if ( this->nLoci_ != this->plaf_.size() ){
        throw LociNumberUnequal( this->plafFileName_ );
    }

    if ( this->nLoci_ != this->altCount_.size() ){
        throw LociNumberUnequal( this->altFileName_ );
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


void DEploidIO::removeFilesWithSameName(){
    strExportProp = this->prefix_ + ".prop";
    strExportLLK = this->prefix_ + ".llk";
    strExportHap = this->prefix_ + ".hap";

    strIbdExportProp = this->prefix_ + ".ibd.prop";
    strIbdExportLLK = this->prefix_ + ".ibd.llk";
    strIbdExportHap = this->prefix_ + ".ibd.hap";
    strIbdExportProbs = this->prefix_ + ".ibd.probs";

    strExportVcf = this->prefix_ + ".vcf";
    if ( compressVcf() ){
        strExportVcf += ".gz";
    }
    strExportLog =  this->prefix_ + ((this->doLsPainting()) ? ".painting":"") + ".log";
    strExportRecombProb = this->prefix_ + ".recomb";

    strExportExtra = this->prefix_ + ".extra";

    if ( this->doLsPainting() == false ){
        if (this->useIBD()){
            remove(strIbdExportProp.c_str());
            remove(strIbdExportLLK.c_str());
            remove(strIbdExportHap.c_str());
        }
        remove(strExportLLK.c_str());
        remove(strExportHap.c_str());
        remove(strExportVcf.c_str());
        remove(strExportProp.c_str());
        remove(strExportExtra.c_str());
        remove(strIbdExportProbs.c_str());
    }

    if (this->doLsPainting() || this->doExportPostProb() ){
        if (this->useIBD()){
            strIbdExportSingleFwdProbPrefix = this->prefix_ + ".ibd.single";
            for ( size_t i = 0; i < this->kStrain_ ; i++ ){
                string tmpStrExportSingleFwdProb = strIbdExportSingleFwdProbPrefix + to_string(i);
                remove(tmpStrExportSingleFwdProb.c_str());
            }
            strIbdExportPairFwdProb = this->prefix_ + ".ibd.pair";
            remove(strIbdExportPairFwdProb.c_str());
        }
        strExportSingleFwdProbPrefix = this->prefix_ + ".single";
        for ( size_t i = 0; i < this->kStrain_ ; i++ ){
            string tmpStrExportSingleFwdProb = strExportSingleFwdProbPrefix + to_string(i);
            remove(tmpStrExportSingleFwdProb.c_str());
        }
        strExportPairFwdProb = this->prefix_ + ".pair";
        remove(strExportPairFwdProb.c_str());
    }
    remove(strExportLog.c_str());
    remove(strExportRecombProb.c_str());

}


void DEploidIO::parse (){
    do {
        if (*argv_i == "-ref") {
            if ( this->useVcf() ){
                throw ( FlagsConflict((*argv_i) , "-vcf") );
            }
            this->readNextStringto ( this->refFileName_ ) ;
        } else if (*argv_i == "-alt") {
            if ( this->useVcf() ){
                throw ( FlagsConflict((*argv_i) , "-vcf") );
            }
            this->readNextStringto ( this->altFileName_ ) ;
        } else if (*argv_i == "-vcf") {
            if ( this->refFileName_.size() > 0 || this->altFileName_.size() > 0 ){
                throw ( FlagsConflict((*argv_i) , "-ref or -alt") );
            }
            this->setUseVcf(true);
            this->readNextStringto ( this->vcfFileName_ ) ;
        } else if (*argv_i == "-vcfOut"){
            this->setDoExportVcf (true);
        } else if (*argv_i == "-plaf") {
            this->readNextStringto ( this->plafFileName_ ) ;
        } else if (*argv_i == "-panel") {
            if ( this->usePanel() == false ){
                throw ( FlagsConflict((*argv_i) , "-noPanel") );
            }
            this->readNextStringto ( this->panelFileName_ ) ;
        } else if (*argv_i == "-noPanel"){
            if ( usePanel() && this->panelFileName_.size() > 0 ){
                throw ( FlagsConflict((*argv_i) , "-panel") );
            }
            if ( doExportPostProb() ){
                throw ( FlagsConflict((*argv_i) , "-exportPostProb") );
            }
            if ( doAllowInbreeding() ){
                throw ( FlagsConflict((*argv_i) , "-inbreeding") );
            }
            this->setUsePanel(false);
            this->setDoExportSwitchMissCopy ( false );
        } else if (*argv_i == "-exclude"){
            this->setExcludeSites( true );
            this->readNextStringto ( this->excludeFileName_ ) ;
        } else if (*argv_i == "-o") {
            this->readNextStringto ( this->prefix_ ) ;
        } else if ( *argv_i == "-p" ) {
            this->precision_ = readNextInput<size_t>() ;
        } else if ( *argv_i == "-k" ) {
            this->setKStrainWasManuallySet(true);
            this->setKstrain(readNextInput<size_t>());

            if ( this->kStrainWasSetByHap() && this->kStrain() != this->initialHap[0].size() ){
                string hint = string(" k = ") + to_string(this->kStrain()) + ", " + this->initialHapFileName_ + " suggests otherwise";
                throw NumOfPropNotMatchNumStrain(hint);
            }

            if ( this->initialPropWasGiven() && this->kStrain_ != initialProp.size() ){
                string hint = string(" k = ") + to_string(kStrain_) + ", flag -initialP suggests otherwise";;
                throw NumOfPropNotMatchNumStrain(hint);
            }

        } else if ( *argv_i == "-nSample" ) {
            this->nMcmcSample_ = readNextInput<size_t>() ;
        } else if ( *argv_i == "-burn" ) {
            this->mcmcBurn_ = readNextInput<double>() ;
            if ( this->mcmcBurn_ < 0 || this->mcmcBurn_ > 1){
                throw ( OutOfRange ("-burn", *argv_i) );
            }
        } else if ( *argv_i == "-miss" ) {
            this->missCopyProb_ = readNextInput<double>() ;
            if ( this->missCopyProb_ < 0 || this->missCopyProb_ > 1){
                throw ( OutOfRange ("-miss", *argv_i) );
            }
        } else if ( *argv_i == "-c" ) {
            this->scalingFactor_ = readNextInput<double>() ;
        } else if ( *argv_i == "-G" ) {
            this->setParameterG(readNextInput<double>());
        } else if ( *argv_i == "-sigma" ) {
            this->setParameterSigma(readNextInput<double>());
        } else if ( *argv_i == "-ibdSigma" ) {
            this->setIBDSigma(readNextInput<double>());
        } else if ( *argv_i == "-recomb" ) {
            this->constRecombProb_ = readNextInput<double>();
            this->useConstRecomb_ = true;
            if ( this->constRecombProb_ < 0 || this->constRecombProb_ > 1){
                throw ( OutOfRange ("-recomb", *argv_i) );
            }
        } else if ( *argv_i == "-printRecomb" ) {
            this->setDoExportRecombProb( true );
        } else if ( *argv_i == "-forbidSame" ) {
            this->setForbidCopyFromSame( true );
        } else if ( *argv_i == "-rate" ) {
            this->mcmcMachineryRate_ = readNextInput<size_t>() ;
        } else if ( *argv_i == "-forbidUpdateProp" ) {
            this->setDoUpdateProp( false );
        } else if ( *argv_i == "-forbidUpdateSingle" ) {
            this->setDoUpdateSingle( false );
        } else if ( *argv_i == "-forbidUpdatePair" ) {
            this->setDoUpdatePair( false );
        } else if ( *argv_i == "-inbreeding" ) {
            if ( this->usePanel() == false ){
                throw ( FlagsConflict((*argv_i) , "-noPanel") );
            }
            this->setDoAllowInbreeding( true );
            this->setDoExportPostProb( true );
        } else if ( *argv_i == "-exportPostProb" ) {
            if ( this->usePanel() == false ){
                throw ( FlagsConflict((*argv_i) , "-noPanel") );
            }
            this->setDoExportPostProb( true );
        } else if ( *argv_i == "-painting" ) {
            if ( this->usePanel() == false ){
                throw ( FlagsConflict((*argv_i) , "-noPanel") );
            }
            if ( this->initialHapWasGiven() == true ){
                throw ( FlagsConflict((*argv_i) , "-initialHap") );
            }
            this->readNextStringto ( this->initialHapFileName_ ) ;
            this->setDoLsPainting( true );
            this->readInitialHaps();
        } else if ( *argv_i == "-ibd" ){
            this->setUseIBD(true);
        } else if ( *argv_i == "-computeLLK" ){
            this->setDoComputeLLK( true );
        } else if ( *argv_i == "-ibdPainting" ){
            this->setDoIbdPainting( true );
        } else if ( *argv_i == "-skipCheckingInitialP" ){
            this->setPleaseCheckInitialP(false);
        } else if ( *argv_i == "-initialP" ){
            this->readInitialProportions();
            this->setInitialPropWasGiven( true );

            // If the k was set manually, check
            if ( this->kStrainWasManuallySet() && this->kStrain_ != initialProp.size() ){
                string hint = string(" k = ") + to_string(kStrain_);
                throw NumOfPropNotMatchNumStrain(hint);
            }

            // If the k was set by initial Hap, check
            if ( this->kStrainWasSetByHap() && this->kStrain() != this->initialProp.size() ){
                string hint = string(" k = ") + to_string(this->kStrain()) + ", " + this->initialHapFileName_ + " suggests otherwise";
                throw NumOfPropNotMatchNumStrain(hint);
            }

            this->setKstrain(this->initialProp.size());
        } else if ( *argv_i == "-initialHap" ){
            if ( this->doLsPainting() == true ){
                throw ( FlagsConflict((*argv_i) , "-painting") );
            }
            this->readNextStringto ( this->initialHapFileName_ ) ;
            this->setInitialHapWasGiven(true);
            this->readInitialHaps();
        } else if ( *argv_i == "-seed"){
            this->set_seed( readNextInput<size_t>() );
            this->setrandomSeedWasGiven( true );
        } else if ( *argv_i == "-z" ){
            this->setCompressVcf(true);
        } else if ( *argv_i == "-h" || *argv_i == "-help"){
            this->set_help(true);
        } else if ( *argv_i == "-v" || *argv_i == "-version"){
            this->setVersion(true);
        } else {
            throw ( UnknowArg((*argv_i)) );
        }
    } while ( ++argv_i != argv_.end());
}


void DEploidIO::checkInput(){

    if ( this->refFileName_.size() == 0 && this->useVcf() == false ){
        throw FileNameMissing ( "Ref count" );}
    if ( this->altFileName_.size() == 0 && this->useVcf() == false ){
        throw FileNameMissing ( "Alt count" );}
    if ( this->plafFileName_.size() == 0 ){
        throw FileNameMissing ( "PLAF" );}
    if ( usePanel() && this->panelFileName_.size() == 0 && !this->doIbdPainting() && !this->doComputeLLK() ){
        throw FileNameMissing ( "Reference panel" );}
    if ( this->initialPropWasGiven() && ( abs(sumOfVec(initialProp) - 1.0) > 0.00001 ) && this->pleaseCheckInitialP() ){
        throw SumOfPropNotOne ( to_string(sumOfVec(initialProp)) );}
    if ( this->initialPropWasGiven() ){
        if ( this->kStrainWasManuallySet() == true ){
        } else {
            // set k strain by proportion length
        }
    }
}


void DEploidIO::readInitialProportions(){
    string tmpFlag = *argv_i;
    ++argv_i;
    if (argv_i == argv_.end() || (*argv_i)[0] == '-' ) {
            throw NotEnoughArg (tmpFlag);
    }

    do {
        try {
            double tmp = convert<double>(*argv_i);
            this->initialProp.push_back(tmp);
        } catch (WrongType e) {
            --argv_i;
            break;
        }
        ++argv_i;
    } while ( argv_i != argv_.end() && (*argv_i)[0] != '-' );
    --argv_i;

    return;
}


void DEploidIO::readNextStringto( string &readto ){
    string tmpFlag = *argv_i;
    ++argv_i;
    if (argv_i == argv_.end() || (*argv_i)[0] == '-' ){
        throw NotEnoughArg(tmpFlag);
    }
    readto = *argv_i;
}


void DEploidIO::printVersion(std::ostream& out){
    out << endl
        << "dEploid " << VERSION
        << endl
        << "Git commit: " << dEploidGitVersion_ << endl;
}

void DEploidIO::printHelp(std::ostream& out){
    out << endl
        << "dEploid " << VERSION
        << endl
        << endl;
    out << "Usage:"
        << endl;
    out << setw(20) << "-h or -help"         << "  --  " << "Help. List the following content."<<endl;
    out << setw(20) << "-v or -version"      << "  --  " << "DEploid version."<<endl;
    out << setw(20) << "-vcf STR"            << "  --  " << "VCF file path."<<endl;
    out << setw(20) << "-ref STR"            << "  --  " << "File path of reference allele count."<<endl;
    out << setw(20) << "-alt STR"            << "  --  " << "File path of alternative allele count."<<endl;
    out << setw(20) << "-plaf STR"           << "  --  " << "File path of population level allele frequencies."<<endl;
    out << setw(20) << "-panel STR"          << "  --  " << "File path of the reference panel."<<endl;
    out << setw(20) << "-exclude STR"        << "  --  " << "File path of sites to be excluded."<<endl;
    out << setw(20) << "-o STR"              << "  --  " << "Specify the file name prefix of the output."<<endl;
    out << setw(20) << "-p INT"              << "  --  " << "Out put precision (default value 8)."<<endl;
    out << setw(20) << "-k INT"              << "  --  " << "Number of strain (default value 5)."<<endl;
    out << setw(20) << "-seed INT"           << "  --  " << "Random seed."<<endl;
    out << setw(20) << "-nSample INT"        << "  --  " << "Number of MCMC samples."<<endl;
    out << setw(20) << "-rate INT"           << "  --  " << "MCMC sample rate."<<endl;
    out << setw(20) << "-noPanel"            << "  --  " << "Use population level allele frequency as prior."<<endl;
    out << setw(20) << "-forbidUpdateProp"   << "  --  " << "Forbid MCMC moves to update proportions."<<endl;
    out << setw(20) << "-forbidUpdateSingle" << "  --  " << "Forbid MCMC moves to update single haplotype."<<endl;
    out << setw(20) << "-forbidUpdatePair"   << "  --  " << "Forbid MCMC moves to update pair haplotypes."<<endl;
    out << setw(20) << "-initialP FLT ..."   << "  --  " << "Initialize proportions."<<endl;
    out << endl;
    out << "Examples:" << endl;
    out << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-CNopanel -noPanel"<< endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf -exclude data/testData/labStrains.test.exclude.txt -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-CNopanelExclude -noPanel"<< endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf -exclude data/testData/labStrains.test.exclude.txt -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-CPanelExclude -panel data/testData/labStrains.test.panel.txt" << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf -exclude data/testData/labStrains.test.exclude.txt -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-CPanelExclude -panel data/testData/labStrains.test.panel.txt -painting PG0390-CPanelExclude.hap" << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-CNopanel -noPanel -k 2 -ibd -nSample 250 -rate 8 -burn 0.67" <<endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-CNopanel -ibdPainting -initialP 0.2 0.8" <<endl;
}


std::ostream& operator<< (std::ostream& stream, const DEploidIO& dEploidIO) {
  for (std::string arg : dEploidIO.argv_) stream << " " << arg;
  return stream;
}


void DEploidIO::readInitialHaps(){
    assert( this->initialHap.size() == 0 );
    InitialHaplotypes initialHapToBeRead;
    initialHapToBeRead.readFromFile(this->initialHapFileName_.c_str());

    assert (this->initialHap.size() == 0 );
    this->initialHap = initialHapToBeRead.content_;

    if ( this->kStrainWasManuallySet() && this->kStrain()!= initialHapToBeRead.truePanelSize() ){
        string hint = string(" k = ") + to_string(this->kStrain()) + ", " + this->initialHapFileName_ + " suggests otherwise";
        throw NumOfPropNotMatchNumStrain(hint);
    }

    if ( this->kStrainWasSetByProp() && this->kStrain() != initialHapToBeRead.truePanelSize() ){
        string hint = string(" k = ") + to_string(kStrain_) + ", flag -initialP suggests otherwise";;
        throw NumOfPropNotMatchNumStrain(hint);
    }

    this->setKstrain(initialHapToBeRead.truePanelSize());
    this->setKStrainWasSetByHap(true);
}


vector <double> DEploidIO::computeExpectedWsafFromInitialHap(){
    // Make this a separate function
    // calculate expected wsaf
    vector <double> expectedWsaf (this->initialHap.size(), 0.0);
    for ( size_t i = 0; i < this->initialHap.size(); i++ ){
        assert( kStrain_ == this->initialHap[i].size() );
        for ( size_t k = 0; k < this->kStrain_; k++){
            expectedWsaf[i] += this->initialHap[i][k] * finalProp[k];
        }
        assert ( expectedWsaf[i] >= 0 );
        //assert ( expectedWsaf[i] <= 1.0 );
    }
    return expectedWsaf;
}


void DEploidIO::computeLLKfromInitialHap(){
    for ( auto const& value: this->initialProp ){
        this->finalProp.push_back(value);
    }

    vector <double> expectedWsaf = computeExpectedWsafFromInitialHap();
    if (expectedWsaf.size() != this->refCount_.size()){
        throw LociNumberUnequal("Hap length differs from data!");
    }
    vector <double> llk = calcLLKs ( this->refCount_, this->altCount_, expectedWsaf, 0, expectedWsaf.size(), this->scalingFactor());
    this->llkFromInitialHap_ = sumOfVec(llk);
}


void DEploidIO::chromPainting(){
    dout << "Painting haplotypes in" << this->initialHapFileName_ <<endl;

    if ( this->initialPropWasGiven() == false ){
        clog << "Initial proportion was not specified. Set even proportions" << endl;
        double evenProp = 1.0 / (double)this->kStrain();
        for ( size_t i = 0; i < this->kStrain(); i++){
            this->initialProp.push_back(evenProp);
        }
    }

    for ( auto const& value: this->initialProp ){
        this->finalProp.push_back(value);
    }

    // Painting posterior probabilities

    // Export the p'
    // Make this a separate class
    //vector < vector <double> > hap;
    //for ( size_t siteI = 0; siteI < decovolutedStrainsToBeRead.content_.size(); siteI++ ){
        //vector <double> tmpHap;
        //for ( size_t tmpk = 0; tmpk < this->kStrain_; tmpk++ ){
            //tmpHap.push_back(decovolutedStrainsToBeRead.content_[siteI][tmpk]);
        //}
        //hap.push_back(tmpHap);
    //}

    //vector < vector <double>> hap = decovolutedStrainsToBeRead.content_;

    vector <double> expectedWsaf = computeExpectedWsafFromInitialHap();

    MersenneTwister tmpRg(this->randomSeed());

    if ( this->doAllowInbreeding() == true ){
        this->panel->initializeUpdatePanel(this->panel->truePanelSize()+kStrain_-1);
    }

    for ( size_t tmpk = 0; tmpk < this->kStrain_; tmpk++ ){
        if ( this->doAllowInbreeding() == true ){
            this->panel->updatePanelWithHaps(this->panel->truePanelSize()+kStrain_-1, tmpk, this->initialHap);
        }

        for ( size_t chromi = 0 ; chromi < this->indexOfChromStarts_.size(); chromi++ ){
            size_t start = this->indexOfChromStarts_[chromi];
            size_t length = this->position_[chromi].size();
            dout << "Painting Chrom "<< chromi << " from site "<< start << " to " << start+length << endl;

            UpdateSingleHap updatingSingle( this->refCount_,
                                      this->altCount_,
                                      this->plaf_,
                                      expectedWsaf,
                                      this->finalProp, this->initialHap, &tmpRg,
                                      start, length,
                                      this->panel, this->missCopyProb_, this->scalingFactor(),
                                      tmpk);

            if ( this->doAllowInbreeding() == true ){
                updatingSingle.setPanelSize(this->panel->inbreedingPanelSize());
            }
            updatingSingle.painting( refCount_, altCount_, expectedWsaf, this->finalProp, this->initialHap);
            //this->writeLastSingleFwdProb( updatingSingle.fwdProbs_, chromi, tmpk, false ); // false as not using ibd
            this->writeLastSingleFwdProb( updatingSingle.fwdBwdProbs_, chromi, tmpk, false ); // false as not using ibd
        }
    }
}


void DEploidIO::readPanel(){
    if ( this->usePanel() == false ){
        return;
    }
    if ( this->doIbdPainting() | this->doComputeLLK() ){
        return;
    }

    panel = new Panel();
    panel->readFromFile(this->panelFileName_.c_str());
    if ( this->excludeSites() ){
        panel->findAndKeepMarkers( this->excludedMarkers );
    }

    panel->computeRecombProbs( this->averageCentimorganDistance(), this->parameterG(), this->useConstRecomb(), this->constRecombProb(), this->forbidCopyFromSame() );
    panel->checkForExceptions( this->nLoci(), this->panelFileName_ );
}


DEploidIO::DEploidIO(const DEploidIO &currentDEploidIO){
    // This is not working! to be improved
    //cout << this->refCount_.size() << endl;
    this->refCount_ = currentDEploidIO.refCount_;
    //cout << this->refCount_.size() << endl;
}


void DEploidIO::getIBDprobsIntegrated(vector < vector <double> > &prob){
    if (prob.size() !=  this->nLoci()){
        throw InvalidInput("Invalid probabilities! Check size!");
    }

    assert(this->ibdProbsIntegrated.size() == 0);

    for (size_t i = 0; i < prob[0].size(); i++){
        this->ibdProbsIntegrated.push_back(0.0);
    }

    for ( size_t siteIndex = 0; siteIndex < this->nLoci(); siteIndex++ ){
        for (size_t i = 0; i < prob[siteIndex].size(); i++){
            this->ibdProbsIntegrated[i] += prob[siteIndex][i];
        }
    }
    normalizeBySum(this->ibdProbsIntegrated);
}


void DEploidIO::computeEffectiveKstrain(vector <double> proportion){
    double tmpSumSq = 0.0;
    for (double p : proportion){
        tmpSumSq += p * p;
    }
    this->effectiveKstrain_ = 1.0 / tmpSumSq;
}


void DEploidIO::computeInferredKstrain(vector <double> proportion){
    this->inferredKstrain_ = 0;
    for (double p : proportion){
        if ( p > 0.01 ){
            this->inferredKstrain_ += 1;
        }
    }
}


void DEploidIO::computeAdjustedEffectiveKstrain(){
    this->adjustedEffectiveKstrain_ = this->effectiveKstrain_;
    if ( (this->inferredKstrain_ == 2) & (ibdProbsIntegrated.size() == 2)){
        if ( this->ibdProbsIntegrated[1] > 0.95 ){
            this->adjustedEffectiveKstrain_ = 1;
        }
    }
}



void DEploidIO::paintIBD(){
    vector <double> goodProp;
    vector <size_t> goodStrainIdx;

    if ( this->doIbdPainting() ){
        this->finalProp = this->initialProp;
    }

    for ( size_t i = 0; i < this->finalProp.size(); i++){
        if (this->finalProp[i] > 0.01){
            goodProp.push_back(this->finalProp[i]);
            goodStrainIdx.push_back(i);
        }
    }

    if (goodProp.size() == 1){
        return;
    }

    DEploidIO tmpDEploidIO; // (*this);
    tmpDEploidIO.setKstrain(goodProp.size());
    tmpDEploidIO.setInitialPropWasGiven(true);
    tmpDEploidIO.initialProp = goodProp;
    tmpDEploidIO.finalProp = goodProp;
    tmpDEploidIO.refCount_ = this->refCount_;
    tmpDEploidIO.altCount_ = this->altCount_;
    tmpDEploidIO.plaf_ = this->plaf_;
    tmpDEploidIO.nLoci_= this->nLoci();
    tmpDEploidIO.position_ = this->position_;
    tmpDEploidIO.chrom_ = this->chrom_;
    //tmpDEploidIO.useConstRecomb_ = true;
    //tmpDEploidIO.constRecombProb_ = 0.000001;

    //tmpDEploidIO.writeLog (&std::cout);

    MersenneTwister tmpRg(this->randomSeed());
    IBDpath tmpIBDpath;
    tmpIBDpath.init(tmpDEploidIO, &tmpRg);
    tmpIBDpath.buildPathProbabilityForPainting(goodProp);
    this->ibdLLK_ = tmpIBDpath.bestPath(goodProp);
    this->ibdProbsHeader = tmpIBDpath.getIBDprobsHeader();
    this->getIBDprobsIntegrated(tmpIBDpath.fwdbwd);
    this->writeIBDpostProb(tmpIBDpath.fwdbwd, this->ibdProbsHeader);
}

