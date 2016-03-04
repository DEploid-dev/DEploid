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

#include <sstream>      // std::stringstream
#include <stdlib.h>     /* strtol, strtod */
#include <fstream>
#include <stdexcept> // std::invalid_argument
#include <vector>
#include <iostream> // std::cout
#include <sstream>      // std::stringstream
#include "global.h"
#include "exceptions.hpp"

#ifndef PARAM
#define PARAM

using namespace std;

class McmcSample;

class PfDeconvIO{
 friend class McmcMachinery;
 friend class TestIO;
  public:
    PfDeconvIO(int argc, char *argv[]);
    ~PfDeconvIO ();

    void init();
    void printHelp();
    bool help() const { return help_; }
    bool usePanel() const { return usePanel_; }
    string panelFileName_;

    // log
    void write (McmcSample * mcmcSample);
    void writeLLK (McmcSample * mcmcSample);
    void writeProp (McmcSample * mcmcSample);
    void writeHap (McmcSample * mcmcSample);
    void writeLog (McmcSample * mcmcSample, ostream * writeTo );

  private:
    // Read in input
    string plafFileName_;
    string refFileName_;
    string altFileName_;
    string prefix_;
    size_t random_seed_;
    bool seed_set_;
    size_t kStrain_;
    size_t precision_;
    size_t nMcmcSample_;
    size_t mcmcMachineryRate_;

    vector <double> plaf_;
    vector <double> refCount_;
    vector <double> altCount_;
    size_t nLoci_;
    bool help_;
    bool usePanel_;

    std::vector<std::string> argv_;
    std::vector<std::string>::iterator argv_i;


    // Output stream
    string pfDeconvVersion_;
    string compileTime_;
    string strExportLLK;
    string strExportHap;
    string strExportProp;
    string strExportLog;
    ofstream ofstreamExportLLK;
    ofstream ofstreamExportHap;
    ofstream ofstreamExportProp;
    ofstream ofstreamExportLog;
    void removeFilesWithSameName();


    // Methods
    void set_panel(const bool usePanel) { this->usePanel_ = usePanel; }
    void set_help(const bool help) { this->help_ = help; }
    void set_seed(const size_t seed){ this->random_seed_ = seed; }

    void parse ();
    void checkInput();
    void finalize();

    void readNextStringto( string &readto );

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

};

#endif
