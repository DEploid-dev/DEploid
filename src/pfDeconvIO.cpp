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
#include <cassert>       // assert
#include <iomanip>      // std::setw

PfDeconvIO::PfDeconvIO( const char plafFileName[],
              const char refFileName[],
              const char altFileName[],
              size_t kStrain ){
    plafFileName_ = string(plafFileName);
    altFileName_ = string(altFileName);
    refFileName_ = string(refFileName);
    (void)readFileLines( refFileName_.c_str(), this->refCount_);
    (void)readFileLines( altFileName_.c_str(), this->altCount_);
    (void)readFileLines( plafFileName_.c_str(), this->plaf_);
    this->nLoci_ = refCount_.size();
    assert( this->nLoci_ == this->plaf_.size() );
    assert( this->altCount_.size() == this->nLoci_ );

    this->kStrain_ = kStrain;
}


PfDeconvIO::~PfDeconvIO(){}


PfDeconvIO::PfDeconvIO(int argc, char *argv[]) {
    argv_ = std::vector<std::string>(argv + 1, argv + argc);
    this->init();
    if ( argv_.size() == 0 ) {
        this->set_help(true);
        return;
    }

    this->parse();
    this->checkInput();
    this->finalize();
}


void PfDeconvIO::init() {
    this->seed_set_ = false;
    this->set_seed( 0 );
    this->set_help(false);
    this->precision_ = 8;
    this->prefix_ = "pf3k-pfDeconv";
    this->kStrain_ = 5;
    this->argv_i = argv_.begin();
}


void PfDeconvIO::finalize(){
    (void)readFileLines( refFileName_.c_str(), this->refCount_);
    (void)readFileLines( altFileName_.c_str(), this->altCount_);
    (void)readFileLines( plafFileName_.c_str(), this->plaf_);
    this->nLoci_ = refCount_.size();
    assert( this->nLoci_ == this->plaf_.size() );
    assert( this->altCount_.size() == this->nLoci_ );
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
            this->readNextStringto ( this->panelFileName_ ) ;
        } else if (*argv_i == "-o") {
            this->readNextStringto ( this->prefix_ ) ;
        } else if ( *argv_i == "-p" ) {
            this->precision_ = readNextInput<size_t>() ;
        } else if ( *argv_i == "-k" ) {
            this->kStrain_ = readNextInput<size_t>() ;
        } else if (*argv_i == "-seed"){
            this->random_seed_ = readNextInput<size_t>() ;
            this->seed_set_ = true;
        } else if (*argv_i == "-h" || *argv_i == "-help"){
            this->set_help(true);
        } else {
            throw ( UnknowArg((*argv_i)) );
        }
    } while ( ++argv_i != argv_.end());
}


void PfDeconvIO::checkInput(){
    if ( this->refFileName_.size() == 0 )
        throw FileNameMissing ( "Ref count" );
    if ( this->altFileName_.size() == 0 )
        throw FileNameMissing ( "Alt count" );
    if ( this->plafFileName_.size() == 0 )
        throw FileNameMissing ( "PLAF" );
    if ( this->panelFileName_.size() == 0 )
        throw FileNameMissing ( "Reference panel" );
}


void PfDeconvIO::readNextStringto( string &readto ){
    string tmpFlag = *argv_i;
    ++argv_i;
    if (argv_i == argv_.end() || (*argv_i)[0] == '-' )
        throw NotEnoughArg(tmpFlag);
    readto = *argv_i;
}


void PfDeconvIO::readFileLines(const char inchar[], vector <double> & out_vec){
    ifstream in_file( inchar );
    string tmp_line;
    if ( in_file.good() ){
        getline ( in_file, tmp_line ); // skip the first line, which is the header
        getline ( in_file, tmp_line );
        while ( tmp_line.size() > 0 ){
            size_t field_start = 0;
            size_t field_end = 0;
            size_t field_index = 0;

            while ( field_end < tmp_line.size() ){
                field_end = min ( tmp_line.find('\t',field_start), tmp_line.find('\n', field_start) );
                string tmp_str = tmp_line.substr( field_start, field_end - field_start );

                if ( field_index == 2 ){
                    out_vec.push_back( strtod(tmp_str.c_str(), NULL) );
                }
                field_start = field_end+1;
                field_index++;
          }
          getline ( in_file, tmp_line );
      }
    } else {
        throw InvalidInputFile( string (inchar) );

    }
    in_file.close();
}


void PfDeconvIO::printHelp(){
    cout << endl
         << "pfDeconv " << VERSION
         << endl
         << endl;
    cout << "Usage:"
         << endl;
    cout << setw(20) << "-h or -help"         << "  --  " << "Help. List the following content."<<endl;
    cout << setw(20) << "-ref STR"            << "  --  " << "Path of reference allele count file."<<endl;
    cout << setw(20) << "-alt STR"            << "  --  " << "Path of alternative allele count file."<<endl;
    cout << setw(20) << "-plaf STR"           << "  --  " << "Path of population level allele frequency file."<<endl;
    cout << setw(20) << "-panel STR"          << "  --  " << "Path of reference panel."<<endl;
    cout << setw(20) << "-o STR"              << "  --  " << "Specify the file name prefix of the output."<<endl;
    cout << setw(20) << "-p INT"              << "  --  " << "Out put precision (default value 8)."<<endl;
    cout << setw(20) << "-k INT"              << "  --  " << "Number of strain (default value 5)."<<endl;
    cout << setw(20) << "-seed INT"           << "  --  " << "Random seed."<<endl
         << endl;
    cout << "Examples:"
         << endl
         << endl;
    cout << "pfDeconv -ref tests/PG0390_first100ref.txt -alt tests/PG0390_first100alt.txt -plaf tests/labStrains_first100_PLAF.txt -panel tests/lab_first100_Panel.txt -o tmp1" << endl;
}
