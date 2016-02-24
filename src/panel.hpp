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

#include <vector>

#ifndef PANEL
#define PANEL

using namespace std;

class Panel{
 friend class TestPanel;
 friend class UpdateSingleHap;
 friend class UpdatePairHap;
 friend class UpdateHap;
  private:
    // Members
    // content is a matrix of n.loci by n.strains, i.e. content length is n.loci
    vector < vector < double > > content_;
    vector < double > recombProbs_;
    size_t nLoci_;
    size_t nPanel_;
    vector <string> chrom_;
    int chromInex_;
    vector < vector < double> > position_;
    vector < double > tmpPosition_;

    // Methods
    void extractChrom( string & tmp_str );
    void extractPOS ( string & tmp_str );
    void computeRecombProbs( double averageCentimorganDistance = 15000.0, double Ne = 10.0);

  public:
    Panel(const char inchar[]);
    ~Panel(){};

    // Methods
    void print();
};

#endif
