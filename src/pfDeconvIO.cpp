/*
 * pfDeconv is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016, Sha (Joe) Zhu, Jacob Almagro and Prof. Gil McVean
 *
 * This file is part of pfDeconv.
 *
 * scrm is free software: you can redistribute it and/or modify
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

#include "pfDeconvIO.hpp"
#include "utility.hpp"  // normailize by sum
#include <cassert>       // assert
#include <iomanip>      // std::setw

PfDeconvIO::PfDeconvIO(){
    this->init();
}


PfDeconvIO::~PfDeconvIO(){
    if ( this->excludedMarkers != NULL ){
        delete this->excludedMarkers;
    }
}


void PfDeconvIO::core(int argc, char *argv[]) {
    argv_ = std::vector<std::string>(argv + 1, argv + argc);
    this->argv_i = argv_.begin();

    if ( argv_.size() == 0 ) {
        this->set_help(true);
        return;
    }

    this->reInit(); // Reset to defalt values before parsing
    this->parse();

    this->checkInput();
    this->finalize();
}


void PfDeconvIO::init() {
    this->randomSeedWasSet_ = false;
    this->initialPropWasGiven_ = false;
    this->initialProp.clear();
    this->exclude_sites_ = false;
    this->excludedMarkers = NULL;
    this->set_seed( 0 );
    this->set_help(false);
    this->set_panel(true);
    this->precision_ = 8;
    this->prefix_ = "pf3k-pfDeconv";
    this->kStrain_ = 5;
    this->nMcmcSample_ = 800;
    this->setDoUpdateProp ( true );
    this->setDoUpdatePair ( true );
    this->setDoUpdateSingle ( true );
    this->setDoExportPostProb( false );
    this->mcmcBurn_ = 0.5;
    this->mcmcMachineryRate_ = 5;
    this->missCopyProb_ = 0.01;
    this->useConstRecomb_ = false;
    this->forbidCopyFromSame_ = false;
    this->constRecombProb_ = 1.0;
    this->averageCentimorganDistance_ = 15000.0;
    this->Ne_ = 10.0;

    #ifdef COMPILEDATE
        compileTime_ = COMPILEDATE;
    #else
        compileTime_ = "";
    #endif

    #ifdef PFDECONVVERSION
        pfDeconvVersion_ = PFDECONVVERSION;
    #else
        pfDeconvVersion_ = "";
    #endif
}


void PfDeconvIO::reInit() {
    this->init();
    this->refFileName_.clear();
    this->altFileName_.clear();
    this->plafFileName_.clear();
    this->panelFileName_.clear();
    this->excludeFileName_.clear();
}


void PfDeconvIO::finalize(){
    if ( !this->randomSeedWasSet_ ){
        this->set_seed( (unsigned)(time(0)) );
    }

    if ( this->exclude_sites_ ){
        excludedMarkers = new ExcludeMarker();
        excludedMarkers->readFromFile(excludeFileName_.c_str());
    }

    InputMarker ref;
    ref.readFromFile(refFileName_.c_str());
    if ( this->exclude_sites_ ){
        ref.removeMarkers( excludedMarkers );
    }
    this->refCount_ = ref.info_;

    InputMarker alt;
    alt.readFromFile(altFileName_.c_str());
    if ( this->exclude_sites_ ){
        alt.removeMarkers( excludedMarkers );
    }
    this->altCount_ = alt.info_;

    InputMarker plaf;
    plaf.readFromFile(plafFileName_.c_str());
    if ( this->exclude_sites_ ){
        plaf.removeMarkers( excludedMarkers );
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
}


void PfDeconvIO::removeFilesWithSameName(){
    strExportLLK = this->prefix_ + ".llk";
    strExportHap = this->prefix_ + ".hap";
    strExportProp = this->prefix_ + ".prop";
    strExportLog =  this->prefix_ + ".log";
    strExportRecombProb = this->prefix_ + ".recomb";
    strExportSingleFwdProb0 = this->prefix_ + ".single1";
    strExportSingleFwdProb1 = this->prefix_ + ".single2";
    strExportPairFwdProb = this->prefix_ + ".pair";

    remove(strExportLLK.c_str());
    remove(strExportHap.c_str());
    remove(strExportProp.c_str());
    remove(strExportLog.c_str());
    remove(strExportRecombProb.c_str());
    remove(strExportSingleFwdProb0.c_str());
    remove(strExportSingleFwdProb1.c_str());
    remove(strExportPairFwdProb.c_str());
}


void PfDeconvIO::parse (){
    do {
        if (*argv_i == "-ref") {
            this->readNextStringto ( this->refFileName_ ) ;
        } else if (*argv_i == "-alt") {
            this->readNextStringto ( this->altFileName_ ) ;
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
            this->set_panel(false);
        } else if (*argv_i == "-exclude"){
            this->exclude_sites_ = true;
            this->readNextStringto ( this->excludeFileName_ ) ;
        } else if (*argv_i == "-o") {
            this->readNextStringto ( this->prefix_ ) ;
        } else if ( *argv_i == "-p" ) {
            this->precision_ = readNextInput<size_t>() ;
        } else if ( *argv_i == "-k" ) {
            this->kStrain_ = readNextInput<size_t>() ;
        } else if ( *argv_i == "-nSample" ) {
            this->nMcmcSample_ = readNextInput<size_t>() ;
        } else if ( *argv_i == "-burn" ) {
            this->mcmcBurn_ = readNextInput<double>() ;
            if ( this->mcmcBurn_ < 0 || this->mcmcBurn_ > 1){
                throw ("out of range");
            }
        } else if ( *argv_i == "-miss" ) {
            this->missCopyProb_ = readNextInput<double>() ;
            if ( this->missCopyProb_ < 0 || this->missCopyProb_ > 1){
                throw ("out of range");
            }
        } else if ( *argv_i == "-recomb" ) {
            this->constRecombProb_ = readNextInput<double>();
            this->useConstRecomb_ = true;
            if ( this->constRecombProb_ < 0 || this->constRecombProb_ > 1){
                throw ("out of range");
            }
        } else if ( *argv_i == "-forbidSame" ) {
            this->forbidCopyFromSame_ = true;
        } else if ( *argv_i == "-rate" ) {
            this->mcmcMachineryRate_ = readNextInput<size_t>() ;
        } else if ( *argv_i == "-forbidUpdateProp" ) {
            this->setDoUpdateProp( false );
        } else if ( *argv_i == "-forbidUpdateSingle" ) {
            this->setDoUpdateSingle( false );
        } else if ( *argv_i == "-forbidUpdatePair" ) {
            this->setDoUpdatePair( false );
        } else if ( *argv_i == "-exportPostProb" ) {
            this->setDoExportPostProb( true );
        } else if ( *argv_i == "-initialP" ){
            this->readInitialProportions();
            this->initialPropWasGiven_ = true;
        } else if ( *argv_i == "-seed"){
            this->random_seed_ = readNextInput<size_t>() ;
            this->randomSeedWasSet_ = true;
        } else if ( *argv_i == "-h" || *argv_i == "-help"){
            this->set_help(true);
        } else {
            throw ( UnknowArg((*argv_i)) );
        }
    } while ( ++argv_i != argv_.end());
}


void PfDeconvIO::checkInput(){
    if ( this->refFileName_.size() == 0 ){
        throw FileNameMissing ( "Ref count" );}
    if ( this->altFileName_.size() == 0 ){
        throw FileNameMissing ( "Alt count" );}
    if ( this->plafFileName_.size() == 0 ){
        throw FileNameMissing ( "PLAF" );}
    if ( usePanel() && this->panelFileName_.size() == 0 ){
        throw FileNameMissing ( "Reference panel" );}
    if ( exclude_sites_ && this->excludeFileName_.size() == 0 ){
        throw FileNameMissing ( "Exclude sites" );}
    if ( this->initialPropWasGiven() && ( abs(sumOfVec(initialProp) - 1.0) > 0.000001 )){
        throw SumOfPropNotOne ( to_string(sumOfVec(initialProp)) );}
    if ( this->initialPropWasGiven() && kStrain_ != initialProp.size() ){
        throw NumOfPropNotMatchNumStrain(""); }
    if ( kStrain_ != 2 && this->doExportPostProb() ){
        throw onlyExportPostProbWhenTwoStrains( to_string (kStrain_ ));}
}


void PfDeconvIO::readInitialProportions(){
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


void PfDeconvIO::readNextStringto( string &readto ){
    string tmpFlag = *argv_i;
    ++argv_i;
    if (argv_i == argv_.end() || (*argv_i)[0] == '-' ){
        throw NotEnoughArg(tmpFlag);
    }
    readto = *argv_i;
}


void PfDeconvIO::printHelp(){
    cout << endl
         << "pfDeconv " << VERSION
         << endl
         << endl;
    cout << "Usage:"
         << endl;
    cout << setw(20) << "-h or -help"         << "  --  " << "Help. List the following content."<<endl;
    cout << setw(20) << "-ref STR"            << "  --  " << "File path of reference allele count."<<endl;
    cout << setw(20) << "-alt STR"            << "  --  " << "File path of alternative allele count."<<endl;
    cout << setw(20) << "-plaf STR"           << "  --  " << "File path of population level allele frequencies."<<endl;
    cout << setw(20) << "-panel STR"          << "  --  " << "File path of the reference panel."<<endl;
    cout << setw(20) << "-exclude STR"        << "  --  " << "File path of sites to be excluded."<<endl;
    cout << setw(20) << "-o STR"              << "  --  " << "Specify the file name prefix of the output."<<endl;
    cout << setw(20) << "-p INT"              << "  --  " << "Out put precision (default value 8)."<<endl;
    cout << setw(20) << "-k INT"              << "  --  " << "Number of strain (default value 5)."<<endl;
    cout << setw(20) << "-seed INT"           << "  --  " << "Random seed."<<endl;
    cout << setw(20) << "-nSample INT"        << "  --  " << "Number of MCMC samples."<<endl;
    cout << setw(20) << "-rate INT"           << "  --  " << "MCMC sample rate."<<endl;
    cout << setw(20) << "-noPanel"            << "  --  " << "Use population level allele frequency as prior."<<endl;
    cout << setw(20) << "-forbidUpdateProp"   << "  --  " << "Forbid MCMC moves to update proportions."<<endl;
    cout << setw(20) << "-forbidUpdateSingle" << "  --  " << "Forbid MCMC moves to update single haplotype."<<endl;
    cout << setw(20) << "-forbidUpdatePair"   << "  --  " << "Forbid MCMC moves to update pair haplotypes."<<endl;
    cout << setw(20) << "-initialP FLT ..."   << "  --  " << "Initialize proportions."<<endl;
    cout << endl;
    cout << "Examples:" << endl;
    cout << endl;
    cout << "./pfDeconv -ref labStrains/PG0390_first100ref.txt -alt labStrains/PG0390_first100alt.txt -plaf labStrains/labStrains_first100_PLAF.txt -panel labStrains/lab_first100_Panel.txt -o tmp1" << endl;
    cout << "./pfDeconv -ref labStrains/PG0390_first100ref.txt -alt labStrains/PG0390_first100alt.txt -plaf labStrains/labStrains_first100_PLAF.txt -panel labStrains/lab_first100_Panel.txt -nSample 100 -rate 3" << endl;
    cout << "./pfDeconv -ref tests/testData/refCountForTesting.csv -alt tests/testData/altCountForTesting.csv -plaf tests/testData/plafForTesting.csv -panel tests/testData/panelForTesting.csv -o tmp"<< endl;
    //cout << "./pfDeconv_dbg -ref labStrains/PG0390_first100ref.txt -alt labStrains/PG0390_first100alt.txt -plaf labStrains/labStrains_first100_PLAF.txt -panel labStrains/lab_first100_Panel.txt -nSample 100 -rate 3" << endl;
    //cout << "./pfDeconv_dbg -ref labStrains/PG0390.C_ref.txt -alt labStrains/PG0390.C_alt.txt -plaf labStrains/labStrains_samples_PLAF.txt -panel labStrains/clonalPanel.csv -nSample 500 -rate 5" << endl;
    //cout << "./pfDeconv -ref labStrains/PG0389.C_ref.txt -alt labStrains/PG0389.C_alt.txt -plaf labStrains/labStrains_samples_PLAF.txt -panel labStrains/clonalPanel.csv -exclude labStrains/PG0389.C.exclude.csv" << endl;
}
