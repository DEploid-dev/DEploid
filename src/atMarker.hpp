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


#ifndef ATMARKER
#define ATMARKER

#include <vector>
#include <string>
#include <cassert>


using namespace std;

class AtMarker{
 friend class TestPanel;
 friend class UpdateSingleHap;
 friend class UpdatePairHap;
 friend class UpdateHap;
 friend class Panel;
 friend class PfDeconvIO;
 friend class InputMarker;
  private:
    // Members
    vector <string> chrom_;
    vector < size_t > indexOfChromStarts_;
    vector < vector < double> > position_;
    // info_ only refers to the first column of the content
    vector <double> info_;

    size_t nLoci_;
    size_t nInfoLines_;

    int tmpChromInex_;
    vector < double > tmpPosition_;

    // Methods
    void extractChrom( string & tmp_str );
    void extractPOS ( string & tmp_str );
    void reshapeContentToInfo();
    void getIndexOfChromStarts();

  public:
    // content is a matrix of n.loci by n.strains, i.e. content length is n.loci
    vector < vector < double > > content_;
    AtMarker(const char inchar[]);
    virtual ~AtMarker(){};
};


class ExcludeMarker : public AtMarker {
  // sorting
  public:
    ExcludeMarker(const char inchar[] ):AtMarker(inchar){};
    ~ExcludeMarker(){};
};


class InputMarker : public AtMarker {
 //friend class PfDeconvIO;
  //private:

  public:
    InputMarker(const char inchar[] ):AtMarker(inchar){ };
    void removeMarkers( ExcludeMarker* excludedMarkers );
    virtual ~InputMarker(){};
};

#endif
