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

#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

#ifndef PARAM
#define PARAM

using namespace std;

class Input{
 friend class McmcMachinery;
 friend class TestInput;
  public:
    Input( const char plafFileName[],
           const char refFileName[],
           const char altFileName[],
           size_t kStrain );
    ~Input ();

  private:
    // Read in input
    vector <double> plaf;
    vector <double> refCount;
    vector <double> altCount;
    size_t kStrain_;
    size_t nLoci_;

    void readFileLines(const char inchar[], vector <double> & out_vec);
};

#endif
