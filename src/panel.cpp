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

#include <fstream>
#include <stdexcept>      // std::invalid_argument
#include <iostream>
#include "panel.hpp"

Panel::Panel(const char inchar[]){
    //cout << "Build Panel" << endl;
    ifstream in_file( inchar );
    string tmp_line;
    if ( in_file.good() ){
        getline ( in_file, tmp_line ); // skip the first line, which is the header
        getline ( in_file, tmp_line );
        while ( tmp_line.size() > 0 ){
            size_t field_start = 0;
            size_t field_end = 0;
            size_t field_index = 0;
            vector <double> contentRow;
            while ( field_end < tmp_line.size() ){
                field_end = min ( tmp_line.find(',',field_start), tmp_line.find('\n', field_start) );
                string tmp_str = tmp_line.substr( field_start, field_end - field_start );

                if ( field_index > 1 ){
                    //cout << strtod(tmp_str.c_str(), NULL) <<" ";
                    contentRow.push_back( strtod(tmp_str.c_str(), NULL) );
                }
                field_start = field_end+1;
                field_index++;
            }
            this->content_.push_back(contentRow);
            getline ( in_file, tmp_line );
        }
    } else {
        throw std::invalid_argument("Invalid input file. " + string (inchar) );

    }
    in_file.close();
    //this->print();

    this->nLoci_ = this->content_.size();
    this->nPanel_ = this->content_.back().size();
    //recombProbs_ = vector <double> (this->nLoci_, 0.0);
    recombProbs_ = vector <double> (this->nLoci_, 0.01); // use constant recombination probability
}


void Panel::print(){
    for ( size_t i = 0; i < this->content_.size(); i++){
        for (size_t j = 0; j < this->content_[i].size(); j++){
            cout <<this->content_[i][j]<<" ";
        }
        cout<<endl;
    }
}
