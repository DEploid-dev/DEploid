/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016, Sha (Joe) Zhu, Jacob Almagro and Prof. Gil McVean
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
 *
 */

#include <stdlib.h>     /* strtol, strtod */
#include <fstream>
#include <stdexcept> // std::invalid_argument
#include <vector>
#include <iostream> // std::cout
#include <sstream>      // std::stringstream
#include "global.h"
#include "exceptions.hpp"
#include "panel.hpp"
#include "vcfReader.hpp"


#ifndef PARAM
#define PARAM

using namespace std;

class McmcSample;
class UpdateSingleHap;
class UpdatePairHap;

class DEploidIO{
#ifdef UNITTEST
 friend class TestIO;
 friend class TestMcmcMachinery;
#endif
 friend class McmcMachinery;
  public:
    DEploidIO();
    DEploidIO(const std::string &arg);
    DEploidIO(int argc, char *argv[]);
    ~DEploidIO ();

    void core();
    void printHelp(std::ostream& out);
    bool help() const { return help_; }
    void printVersion(std::ostream& out);
    bool version() const { return version_; }

    bool initialPropWasGiven() const { return initialPropWasGiven_; }

    // Panel related
    bool usePanel() const { return usePanel_; }
    string panelFileName_;

    size_t nLoci() const { return this->nLoci_; }
    double averageCentimorganDistance() const { return this->averageCentimorganDistance_; }
    double Ne() const { return this->Ne_; }
    double constRecombProb() const { return this->constRecombProb_; }
    bool useConstRecomb() const { return this->useConstRecomb_; }

    ExcludeMarker* excludedMarkers;
    bool excludeSites_;
    bool excludeSites() const {return this->excludeSites_; }
    void setExcludeSites(const size_t exclude){ this->excludeSites_ = exclude; }

    bool forbidCopyFromSame() const { return this->forbidCopyFromSame_; }
    void setForbidCopyFromSame(const bool forbid){ this->forbidCopyFromSame_ = forbid; }

    // Log
    void write (McmcSample * mcmcSample, Panel * panel );
    bool randomSeedWasSet() const {return this->randomSeedWasSet_; }

    friend std::ostream& operator<< (std::ostream& stream, const DEploidIO& dEploidIO);
    size_t randomSeed() const { return randomSeed_;}

  private:

    // Read in input
    string plafFileName_;
    string refFileName_;
    string altFileName_;
    string vcfFileName_;
    string excludeFileName_;
    string prefix_;
    size_t randomSeed_;
    bool randomSeedWasSet_;
    void setRandomSeedWasSet(const bool random){ this->randomSeedWasSet_ = random; }


    bool initialPropWasGiven_;
    bool useConstRecomb_;
    bool forbidCopyFromSame_;
    size_t kStrain_;
    size_t precision_;
    size_t nMcmcSample_;
    size_t mcmcMachineryRate_;
    double mcmcBurn_;

    bool doUpdateProp_;
    bool doUpdatePair_;
    bool doUpdateSingle_;
    bool doExportPostProb_;
    bool doExportSwitchMissCopy_;

    vector <double> initialProp;
    vector <string> chrom_;
    vector < size_t > indexOfChromStarts_;
    vector < vector < int > > position_;
    vector <double> plaf_;
    vector <double> refCount_;
    vector <double> altCount_;
    size_t nLoci_;

    // Help related
    bool help_;
    void set_help(const bool help) { this->help_ = help; }
    bool version_;
    void setVersion(const bool version) { this->version_ = version; }

    // Panel related
    bool usePanel_;
    void set_panel(const bool usePanel) { this->usePanel_ = usePanel; }

    // Vcf Related
    VcfReader * vcfReaderPtr_;

    bool useVcf_;
    void setUseVcf(const bool useVcf){ this->useVcf_ = useVcf; }
    bool useVcf() const {return this->useVcf_; }

    bool doExportVcf_;
    void setDoExportVcf( const bool exportVcf ){ this->doExportVcf_ = exportVcf; }
    bool doExportVcf() const { return this->doExportVcf_; }

    bool compressVcf_;
    void setCompressVcf( const bool compress ){ this->compressVcf_ = compress; }
    bool compressVcf() const { return this->compressVcf_; }

    bool doExportRecombProb_;
    void setDoExportRecombProb( const bool exportRecombProb ){ this->doExportRecombProb_ = exportRecombProb; }
    bool doExportRecombProb() const { return this->doExportRecombProb_; }

    // Parameters
    double missCopyProb_;
    double averageCentimorganDistance_;// = 15000.0,
    double Ne_;// = 10.0
    double constRecombProb_;

    std::vector<std::string> argv_;
    std::vector<std::string>::iterator argv_i;


    // Output stream
    string dEploidGitVersion_;
    string compileTime_;
    string strExportLLK;
    string strExportHap;
    string strExportVcf;
    string strExportProp;
    string strExportLog;
    string strExportRecombProb;

    string strExportSingleFwdProbPrefix;
    string strExportPairFwdProb;

    string strExportOneSwitchOne;
    string strExportOneMissCopyOne;
    string strExportTwoSwitchOne;
    string strExportTwoSwitchTwo;
    string strExportTwoMissCopyOne;
    string strExportTwoMissCopyTwo;

    ofstream ofstreamExportTmp;
    ofstream ofstreamExportFwdProb;

    string startingTime_;
    string endTime_;
    void getTime( bool isStartingTime );


    // Methods
    void init();
    void reInit();
    void parse ();
    void checkInput();
    void finalize();
    void readNextStringto( string &readto );
    void readInitialProportions();

    void set_seed(const size_t seed){ this->randomSeed_ = seed; }
    void removeFilesWithSameName();


    template<class T>
    T readNextInput() {
        string tmpFlag = *argv_i;
        ++argv_i;
        if (argv_i == argv_.end() || (*argv_i)[0] == '-' ) {
            throw NotEnoughArg (tmpFlag);}
        return convert<T>(*argv_i);
    }

    template<class T>
    T convert(const std::string &arg) {
        T value;
        std::stringstream ss(arg);
        ss >> value;
        if (ss.fail()) {
            throw WrongType(arg);
        }
        return value;
    }

    // Getters and Setters
    void setDoUpdateProp ( const bool setTo ){ this->doUpdateProp_ = setTo; }
    bool doUpdateProp() const { return this->doUpdateProp_; }

    void setDoUpdateSingle ( const bool setTo ){ this->doUpdateSingle_ = setTo; }
    bool doUpdateSingle() const { return this->doUpdateSingle_; }

    void setDoUpdatePair ( const bool setTo ){ this->doUpdatePair_ = setTo; }
    bool doUpdatePair() const { return this->doUpdatePair_; }

    void setDoExportPostProb ( const bool setTo ){ this->doExportPostProb_ = setTo; }
    bool doExportPostProb() const { return this->doExportPostProb_; }

    void setDoExportSwitchMissCopy ( const bool setTo ){ this->doExportSwitchMissCopy_ = setTo; }
    bool doExportSwitchMissCopy() const { return this->doExportSwitchMissCopy_; }

    // log and export resutls
    void writeRecombProb ( Panel * panel );
    void writeLLK (McmcSample * mcmcSample);
    void writeProp (McmcSample * mcmcSample);
    void writeHap (McmcSample * mcmcSample);
    void writeVcf (McmcSample * mcmcSample);
    void writeLog (McmcSample * mcmcSample, ostream * writeTo );
    void writeLastSingleFwdProb( UpdateSingleHap & updateSingle, size_t chromIndex, size_t strainIndex  );
    void writeLastPairFwdProb( UpdatePairHap & updatePair, size_t chromIndex );
};

#endif
