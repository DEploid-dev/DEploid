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

#include "fastfunc.hpp"

void FastFunc::build_fastlog_double_table(int size) {
  fastlog_double_table_ = std::vector<double>(size+1);
  double prevx = 1.0;
  double prevy = 0.0;
  for (int index=0; index<size+1; index++) {

    // calculate x coordinate at which linear approximation is exactly equal
    // to the true logarithm
    double curx = 1.0 + (index+5.0/6.0) / size;
    if (index == size-1)
      curx = 1.0 + (index+1.0) / size;

    // calculate true logarithm at the point
    double cury = log( curx );

    // calculate the linear approximation at the next join point
    double targetx = 1.0 + (index+1.0)/size;
    double targety = prevy + (targetx-prevx)*(cury-prevy)/(curx-prevx);

    // store previous linear approximation in table, and update prevx/y
    fastlog_double_table_.at(index) = (double)(prevy);
    prevx = targetx;
    prevy = targety;
  }
}
