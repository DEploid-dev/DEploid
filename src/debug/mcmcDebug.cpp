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

#include "mcmc.hpp"
//#include "utility.hpp"
//#include <math.h>       /* ceil */
//#include <random>
//#include "updateHap.hpp"
//#include <stdio.h>


bool McmcMachinery::doutProp(){
    dout << "  Update proportion to: ";

    for ( auto const& value: this->currentProp_ ){
        dout << value << " ";
    }

    dout<<endl;
    return true;
}


bool McmcMachinery::doutLLK(){
    dout << " Current log likelihood = " << sumOfVec( this->currentLLks_ ) << endl;
    return true;
}
