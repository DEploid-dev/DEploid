// GNU General Public License Agreement
// Copyright (C) 2004-2010 CodeCogs, Zyba Ltd, Broadwood, Holford, TA5 1DU, England.
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by CodeCogs. 
// You must retain a copy of this licence in all copies. 
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
// ---------------------------------------------------------------------------------
//! Evaluates the natural logarithm of a number + 1.

#ifndef MATHS_ARITHMETIC_LN_ADD1_H
#define MATHS_ARITHMETIC_LN_ADD1_H

#include <math.h>

namespace Maths
{

namespace Arithmetic
{

//! Evaluates the natural logarithm of a number + 1.

double ln_add1( double x )
{

  if( fabs(x) > 0.375 )
    return log( 1+x );

  const double p1 = -1.29418923021993;
  const double p2 = 0.405303492862024;
  const double p3 = -1.78874546012214;
  const double q1 = -1.62752256355323;
  const double q2 = 0.747811014037616;
  const double q3 = -0.0845104217945565;

  double t = x/(x+2);
  double t2 = t * t;

  return 2.0 * t * ( ( ( p3*t2 + p2 )*t2 + p1 )*t2 + 1.0 ) /
    ( ( ( q3*t2 + q2 )*t2 + q1 )*t2 + 1.0 );
}

}

}

#endif

