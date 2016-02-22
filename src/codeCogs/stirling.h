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
//! Stirling series approximation of the gamma function

#ifndef MATHS_SPECIAL_GAMMA_STIRLING_H
#define MATHS_SPECIAL_GAMMA_STIRLING_H

#include <poly_eval.h>
#include <math.h>
#include <float.h>

namespace Maths
{

namespace Special
{

namespace Gamma
{

//! Stirling series approximation of the gamma function

double stirling(double x)
{
  static double STIR[5] = {
    7.87311395793093628397E-4,    // ~ 163879/209018880
    -2.29549961613378126380E-4,   // ~ -571/2488320
    -2.68132617805781232825E-3,   // ~ -139/51840
    3.47222221605458667310E-3,    // ~ 1/288
    8.33333333333482257126E-2,    // ~ 1/12
  };
  static double SQTPI = 2.50662827463100050242;  // sqrt of 2*pi

  double w = 1.0/x;
  w = 1.0 + w * Maths::Algebra::Polynomial::polyEval( w, STIR, 4 );

  double y = exp(x);

  if( x > 143.01608 )
  { // Avoid overflow in pow()
    double v = pow( x, 0.5 * x - 0.25 );
    y = v * (v / y);
  }
  else
  {
    y = pow( x, x - 0.5 ) / y;
  }
  y = SQTPI * y * w;
  return( y );
}

}

}

}

#endif

