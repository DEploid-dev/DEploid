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
//! Evaluates ln(gamma(b) / gamma(a+b)).

#ifndef MATHS_SPECIAL_GAMMA_LOGGAMMAFRAC_H
#define MATHS_SPECIAL_GAMMA_LOGGAMMAFRAC_H

#include "ln_add1.h"
#include <math.h>

namespace Maths
{

namespace Special
{

namespace Gamma
{

//! Evaluates ln(gamma(b) / gamma(a+b)).

double logGammaFrac( double a, double b )
{
  double c,d,h,x;
  if( a <= b )
  {
    h = a/b;
    c = h / (1.0+h);
    x = 1.0 / (1.0+h);
    d = b + (a - 0.5);
  }
  else
  {
    h = b/a;
    c = 1.0 / (1.0+h);
    x = h / (1.0+h);
    d = a + b - 0.5;
  }
  const double c0 = 0.0833333333333333;
  const double c1 = -0.00277777777760991;
  const double c2 = 0.000793650666825390;
  const double c3 = -0.000595202931351870;
  const double c4 = 0.000837308034031215;
  const double c5 = -0.00165322962780713;

  double x2 = x*x;
  double s3 = 1.0 + x + x2;
  double s5 = 1.0 + x + x2*s3;
  double s7 = 1.0 + x + x2*s5;
  double s9 = 1.0 + x + x2*s7;
  double s11 = 1.0 + x + x2*s9;
  double t = pow( 1.0/b, 2.0 );
  double w = (((( c5*s11*t + c4*s9 )*t + c3*s7)*t + c2*s5)*t + c1*s3)*t + c0;

  return w*c/b - d*Maths::Arithmetic::ln_add1(a/b) - a*( log(b) - 1.0 );
}

}

}

}

#endif

