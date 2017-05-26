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


#ifndef TXTREADER
#define TXTREADER

#include "variantIndex.hpp"
#include "exceptions.hpp"

class TxtReader : public VariantIndex {
#ifdef UNITTEST
 friend class TestPanel;
 friend class TestTxtReader;
 friend class TestInitialHaplotypes;
#endif
 friend class McmcMachinery;
 friend class UpdateSingleHap;
 friend class UpdatePairHap;
 friend class UpdateHap;
 friend class Panel;
 friend class DEploidIO;
  private:
    // Members

    // content is a matrix of n.loci by n.strains, i.e. content length is n.loci
    vector < vector < double > > content_;
    vector < vector < double > > keptContent_;
    // info_ only refers to the first column of the content
    vector <double> info_;

    size_t nInfoLines_;

    int tmpChromInex_;
    vector < int > tmpPosition_;

    // Methods
    void extractChrom( string & tmp_str );
    void extractPOS ( string & tmp_str );
    void reshapeContentToInfo();
    string fileName;

  public: // move the following to private
    TxtReader (){};
    virtual void readFromFile( const char inchar[] ){ this->readFromFileBase( inchar ); };
    void readFromFileBase( const char inchar[] );
    virtual ~TxtReader(){ };
    void removeMarkers ( );
};




class ExcludeMarker : public TxtReader {
  // sorting
  public:
    ExcludeMarker():TxtReader(){};
    ~ExcludeMarker(){};
};


#endif
