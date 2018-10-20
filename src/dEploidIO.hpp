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
 *
 */

#include <vector>
#include <fstream>
#include <stdlib.h>             // strtol, strtod
#include <stdexcept>            // std::invalid_argument
#include <iostream>             // std::cout
#include <sstream>              // std::stringstream
#include "global.hpp"
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
 friend class TestIBDpath;
#endif
 friend class McmcMachinery;
 friend class RMcmcSample;
 friend class IBDpath;
  public:
    DEploidIO();
    DEploidIO(const DEploidIO &currentDEploidIO);
    DEploidIO(const std::string &arg);
    DEploidIO(int argc, char *argv[]);
    ~DEploidIO ();

    void printHelp(std::ostream& out);
    bool help() const { return help_; }
    void printVersion(std::ostream& out);
    bool version() const { return version_; }
    // Painting related
    void chromPainting ();
    bool doLsPainting() const { return this->doLsPainting_; }
    bool doIbdPainting() const { return this->doIbdPainting_; }
    bool doComputeLLK() const { return this->doComputeLLK_; }
    void computeLLKfromInitialHap();
    bool useIBD() const { return this->useIBD_;}
    bool useIbdOnly() const { return this->useIbdOnly_;}
    bool useLasso() const { return this->useLasso_;}
    void paintIBD();
    double ibdLLK_;
    void getIBDprobsIntegrated(vector < vector <double> > &prob);
    // Lasso related
    void dEploidLasso();
    vector <double> finalProp;

    // Log
    void wrapUp();
    bool randomSeedWasSet() const {return this->randomSeedWasGiven_; }

    friend std::ostream& operator<< (std::ostream& stream, const DEploidIO& dEploidIO);

    size_t randomSeed() const { return randomSeed_;}


    vector <double> plaf_;
    vector <double> refCount_;
    vector <double> altCount_;
    Panel *panel;
    vector < size_t > indexOfChromStarts_;
    vector < vector < int > > position_;

    // Lasso Related
    //vector < vector < vector <double> > > lassoPanels;
    vector < Panel* > lassoPanels;
    vector < vector <double> > lassoPlafs;
    vector < vector <double> > lassoRefCount;
    vector < vector <double> > lassoAltCount;

    void writeHap (vector < vector <double> > &hap, bool useIBD = false);
    bool doPrintLassoPanel_;

    // Trimming related
    void ibdTrimming();



    void setInitialHapWasGiven(const bool setTo) { this->initialHapWasGiven_ = setTo; }
    vector < vector <double> > initialHap;
    vector <double> initialProp;
    void setDoUpdateProp ( const bool setTo ) { this->doUpdateProp_ = setTo; }
    void setInitialPropWasGiven(const bool setTo) {this->initialPropWasGiven_ = setTo; }

  private:
    void core();
    double llkFromInitialHap_;

    // Read in input
    string plafFileName_;
    string refFileName_;
    string altFileName_;
    string vcfFileName_;
    string excludeFileName_;
    string initialHapFileName_;
    string prefix_;
    size_t randomSeed_;
    bool randomSeedWasGiven_;
    void setrandomSeedWasGiven(const bool random) { this->randomSeedWasGiven_ = random; }


    bool initialPropWasGiven_;
    bool pleaseCheckInitialP_;
    bool initialHapWasGiven_;
    bool kStrainWasManuallySet_;
    bool kStrainWasSetByHap_;
    bool kStrainWasSetByProp_;
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
    bool doAllowInbreeding_;
    bool doLsPainting_;
    bool doIbdPainting_;
    bool useIBD_;
    bool useIbdOnly_;
    bool useLasso_;

    vector <string> chrom_;
    vector <double> obsWsaf_;
    vector <size_t> wsafGt0At_;
    size_t nLoci_;

    // Help related
    bool help_;
    void set_help(const bool help) { this->help_ = help; }
    bool version_;
    void setVersion(const bool version) { this->version_ = version; }

    // Panel related
    bool usePanel_;
    void setUsePanel(const bool setTo) { this->usePanel_ = setTo; }

    // Vcf Related
    VcfReader * vcfReaderPtr_;

    bool useVcf_;
    void setUseVcf(const bool useVcf) { this->useVcf_ = useVcf; }
    bool useVcf() const {return this->useVcf_; }

    bool doExportVcf_;
    void setDoExportVcf( const bool exportVcf ) { this->doExportVcf_ = exportVcf; }
    bool doExportVcf() const { return this->doExportVcf_; }

    bool compressVcf_;
    void setCompressVcf( const bool compress ) { this->compressVcf_ = compress; }
    bool compressVcf() const { return this->compressVcf_; }

    bool doExportRecombProb_;
    void setDoExportRecombProb( const bool exportRecombProb ) { this->doExportRecombProb_ = exportRecombProb; }
    bool doExportRecombProb() const { return this->doExportRecombProb_; }

    bool doComputeLLK_;
    void setDoComputeLLK( const bool setTo ) { this->doComputeLLK_ = setTo; }



    // Parameters
    double missCopyProb_;
    double averageCentimorganDistance_;// = 15000.0,
    //double Ne_;// = 10.0
    double constRecombProb_;
    double scalingFactor_; // 100.0

    std::vector<std::string> argv_;
    std::vector<std::string>::iterator argv_i;

    // Diagnostics
    double maxLLKs_;
    void setmaxLLKs ( const double setTo ) { this->maxLLKs_ = setTo; }
    double meanThetallks_;
    void setmeanThetallks ( const double setTo ) { this->meanThetallks_ = setTo; }
    double meanllks_;
    void setmeanllks ( const double setTo ) { this->meanllks_ = setTo; }
    double stdvllks_;
    void setstdvllks ( const double setTo ) { this->stdvllks_ = setTo; }
    double dicByTheta_;
    void setdicByTheta ( const double setTo ) { this->dicByTheta_ = setTo; }
    double dicByVar_;
    void setdicByVar ( const double setTo ) { this->dicByVar_ = setTo; }
    double acceptRatio_;
    void setacceptRatio ( const double setTo ) { this->acceptRatio_ = setTo; }

    // Output stream
    string dEploidGitVersion_;
    string lassoGitVersion_;
    string compileTime_;
    string strExportLLK;
    string strExportHap;
    string strExportVcf;
    string strExportProp;
    string strExportLog;
    string strExportRecombProb;

    string strIbdExportProp;
    string strIbdExportLLK;
    string strIbdExportHap;
    string strIbdExportProbs;

    string strExportSingleFwdProbPrefix;
    string strExportPairFwdProb;
    string strIbdExportSingleFwdProbPrefix;
    string strIbdExportPairFwdProb;

    string strExportExtra;

    ofstream ofstreamExportTmp;
    ofstream ofstreamExportFwdProb;

    string startingTime_;
    string endTime_;
    void getTime( bool isStartingTime );

    bool isCopied_;

    // Methods
    void init();
    void reInit();
    void parse ();
    void checkInput();
    void finalize();
    void readNextStringto( string &readto );
    void readInitialProportions();
    void readInitialHaps();

    void set_seed(const size_t seed) { this->randomSeed_ = seed; }
    void removeFilesWithSameName();
    vector <double> computeExpectedWsafFromInitialHap();


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
    void setDoPrintLassoPanel ( const bool setTo ) { this->doPrintLassoPanel_ = setTo; }
    bool doPrintLassoPanel() const { return this->doPrintLassoPanel_; }

    void setIsCopied ( const bool setTo ) { this->isCopied_ = setTo; }
    bool isCopied() const { return this->isCopied_; }

    bool doUpdateProp() const { return this->doUpdateProp_; }

    void setDoUpdateSingle ( const bool setTo ) { this->doUpdateSingle_ = setTo; }
    bool doUpdateSingle() const { return this->doUpdateSingle_; }

    void setDoUpdatePair ( const bool setTo ) { this->doUpdatePair_ = setTo; }
    bool doUpdatePair() const { return this->doUpdatePair_; }

    void setDoExportPostProb ( const bool setTo ) { this->doExportPostProb_ = setTo; }
    bool doExportPostProb() const { return this->doExportPostProb_; }

    void setDoExportSwitchMissCopy ( const bool setTo ) { this->doExportSwitchMissCopy_ = setTo; }
    bool doExportSwitchMissCopy() const { return this->doExportSwitchMissCopy_; }

    void setDoAllowInbreeding ( const bool setTo ) { this->doAllowInbreeding_ = setTo; }
    bool doAllowInbreeding() const { return this->doAllowInbreeding_; }

    void setDoLsPainting ( const bool setTo ) { this->doLsPainting_ = setTo; }
    void setDoIbdPainting ( const bool setTo ) { this->doIbdPainting_ = setTo; }
    void setUseIBD( const bool setTo) { this->useIBD_ = setTo; }
    void setUseIbdOnly(const bool setTo) { this->useIbdOnly_ = setTo;}
    void setUseLasso( const bool setTo) { this->useLasso_ = setTo; }

    bool initialPropWasGiven() const { return initialPropWasGiven_; }

    bool pleaseCheckInitialP() const { return pleaseCheckInitialP_; }
    void setPleaseCheckInitialP(const bool setTo) {this->pleaseCheckInitialP_ = setTo; }

    bool initialHapWasGiven() const { return initialHapWasGiven_; }

    bool randomSeedWasGiven() const {return this->randomSeedWasGiven_; }

    // log and export resutls
    void writeRecombProb ( Panel * panel );
    void writeIBDpostProb(vector < vector <double> > & reshapedProbs, vector <string> header);
    vector <string> ibdProbsHeader;
    vector <double> ibdProbsIntegrated;

    void writeLLK (McmcSample * mcmcSample, bool useIBD = false);
    void writeProp (McmcSample * mcmcSample, bool useIBD = false);
    void writeVcf (McmcSample * mcmcSample);
    void writeLastSingleFwdProb( vector < vector <double> >& probabilities, size_t chromIndex, size_t strainIndex, bool useIBD );
    void writeLastPairFwdProb( UpdatePairHap & updatePair, size_t chromIndex );
    void writeLog (ostream * writeTo );
    void writeEventCount();

    vector <double> IBDpathChangeAt;
    vector <double> finalIBDpathChangeAt;

    vector <double> siteOfTwoSwitchOne;
    vector <double> siteOfTwoMissCopyOne;
    vector <double> siteOfTwoSwitchTwo;
    vector <double> siteOfTwoMissCopyTwo;
    vector <double> siteOfOneSwitchOne;
    vector <double> siteOfOneMissCopyOne;

    vector <double> finalSiteOfTwoSwitchOne;
    vector <double> finalSiteOfTwoMissCopyOne;
    vector <double> finalSiteOfTwoSwitchTwo;
    vector <double> finalSiteOfTwoMissCopyTwo;
    vector <double> finalSiteOfOneSwitchOne;
    vector <double> finalSiteOfOneMissCopyOne;


    void writeMcmcRelated (McmcSample * mcmcSample, bool useIBD = false);
    void readPanel();

    // Panel related
    bool usePanel() const { return usePanel_; }
    string panelFileName_;
    double parameterG_;
    void setParameterG ( const double setTo ) { this->parameterG_ = setTo; }
    double parameterG() const { return this->parameterG_; }
    double parameterSigma_;
    void setParameterSigma ( const double setTo ) { this->parameterSigma_ = setTo; }
    double parameterSigma() const { return this->parameterSigma_; }
    double ibdSigma_;
    void setIBDSigma ( const double setTo ) { this->ibdSigma_ = setTo; }
    double ibdSigma() const {return this->ibdSigma_;}

    void setNLoci ( const size_t setTo ) { this->nLoci_ = setTo;}
    size_t nLoci() const { return this->nLoci_; }
    void setKstrain ( const size_t setTo ) { this->kStrain_ = setTo;}
    size_t kStrain() const { return this->kStrain_;}
    void setKStrainWasManuallySet ( const size_t setTo ) { this->kStrainWasManuallySet_ = setTo; }
    bool kStrainWasSetByHap() const { return this->kStrainWasSetByHap_; }
    void setKStrainWasSetByHap ( const size_t setTo ) { this->kStrainWasSetByHap_ = setTo; }
    bool kStrainWasManuallySet() const { return this->kStrainWasManuallySet_; }
    void setKStrainWasSetByProp ( const size_t setTo ) { this->kStrainWasSetByProp_ = setTo; }
    bool kStrainWasSetByProp() const { return this->kStrainWasSetByProp_; }
    size_t nMcmcSample() const { return this->nMcmcSample_; }
    double averageCentimorganDistance() const { return this->averageCentimorganDistance_; }
    //double Ne() const { return this->Ne_; }
    double scalingFactor() const {return this->scalingFactor_; }
    void setScalingFactor ( const double setTo ) { this->scalingFactor_ = setTo; }
    double constRecombProb() const { return this->constRecombProb_; }
    bool useConstRecomb() const { return this->useConstRecomb_; }

    ExcludeMarker* excludedMarkers;
    bool excludeSites_;
    bool excludeSites() const {return this->excludeSites_; }
    void setExcludeSites(const size_t exclude) { this->excludeSites_ = exclude; }

    bool forbidCopyFromSame() const { return this->forbidCopyFromSame_; }
    void setForbidCopyFromSame(const bool forbid) { this->forbidCopyFromSame_ = forbid; }

    double effectiveKstrain_ ;
    void computeEffectiveKstrain(vector <double> proportion);
    int inferredKstrain_;
    void computeInferredKstrain(vector <double> proportion);
    double adjustedEffectiveKstrain_;
    void computeAdjustedEffectiveKstrain();

    // Lasso related
    vector <double> lassoComputeObsWsaf(size_t segmentStartIndex, size_t nLoci);
    vector < vector <double> > lassoSubsetPanel(size_t segmentStartIndex, size_t nLoci);
    void writePanel(Panel *panel, size_t chromi, vector <string> hdr);


};

#endif
