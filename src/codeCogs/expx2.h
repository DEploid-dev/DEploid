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
//! Computes the exponential of a squared argument.

#ifndef MATHS_ARITHMETIC_EXPX2_H
#define MATHS_ARITHMETIC_EXPX2_H

#include <math.h>
#include <float.h>

#define MAXLOG 7.08396418532264106224E2
#define M 128.0
#define MINV .0078125

namespace Maths
{

namespace Arithmetic
{

//! Computes the exponential of a squared argument.

double expx2( double x, int sign=1)
{
  x = fabs(x);
  if (sign < 0)
    x = -x;

/* Represent x as an exact multiple of M plus a residual.
     M is a power of 2 chosen so that exp(m * m) does not overflow
     or underflow and so that |x - m| is small.  */
  double m = MINV * floor(M * x + 0.5);
  double f = x - m;

  // x^2 = m^2 + 2mf + f^2
  double u = m * m;
  double u1 = 2*m*f  +  f*f;
  if (sign < 0)
   {
      u = -u;
      u1 = -u1;
    }

  if ((u+u1) > MAXLOG)
    return DBL_MAX;        // from float.h

  // u is exact, u1 is small.
  u = exp(u) * exp(u1);
  return(u);
}

}

}

#undef MAXLOG
#undef M
#undef MINV

#endif

