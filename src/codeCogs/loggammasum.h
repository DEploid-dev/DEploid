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
//! Returns the natural logarithm of the gamma function of the sum of two numbers.

#ifndef MATHS_SPECIAL_GAMMA_LOGGAMMASUM_H
#define MATHS_SPECIAL_GAMMA_LOGGAMMASUM_H

#include <assert.h>
#include <log_gamma.h>
#include <math.h>

namespace Maths
{

namespace Special
{

namespace Gamma
{

//! Returns the natural logarithm of the gamma function of the sum of two numbers.

inline double logGammaSum( double a, double b )
{
  assert( a>=1 && a<=2 );
  assert( b>=1 && b<=2 );
  return Maths::Special::Gamma::log_gamma(a+b);
}

}

}

}

#endif

