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
//! Returns the natural logarithm of the complete beta function.

#ifndef MATHS_SPECIAL_GAMMA_LOGBETA_H
#define MATHS_SPECIAL_GAMMA_LOGBETA_H

#include <math.h>
#include <ln_add1.h>
#include <log_gamma.h>
#include <loggammasum.h>
#include <loggammafrac.h>
#include <asympt_expn.h>

namespace Maths
{

namespace Special
{

namespace Gamma
{

//! Returns the natural logarithm of the complete beta function.

double logBeta( double a, double b )
{
  assert(a>=0 && b>=0);
  const double e = 0.918938533204673;

  if( a > b )
  {
    double t=a;
    a=b;
    b=t;
  }

  if( a >= 8 )
  {
    double h = a/b;
    return e + Maths::Algebra::Series::bcorr(a,b) + ( a-0.5 )*log( h/(1.0+h) )
      - 0.5*log(b) - b * Maths::Arithmetic::ln_add1(h);
  }
  if( a >= 1 )
  {
    if( a > 2 )
    {
      if( b > 1000 )
      {
        double n = a - 1.0;
        double w = 1.0;
        for( int i=1; i<=n; i++ )
        {
          a -= 1.0;
          w *= a / (1.0 + a/b);
        }
        return log(w) - n * log(b) + log_gamma(a) + logGammaFrac(a,b);
      }

      int n = int(a - 1);
      double w = 1.0;
      for( int i=1; i<=n; i++ )
      {
        a -= 1.0;
        double h = a/b;
        w *= h/(1.0+h);
      }

      w = log(w);
      if( b < 8 )
      {
        n = int(b - 1);
        double z = 1.0;
        for( int i=1; i<=n; i++ )
        {
          b -= 1.0;
          z *= b/(a+b);
        }
        return w + log(z) + log_gamma(a) + log_gamma(b) - logGammaSum(a,b);
      }
      return w + log_gamma(a) + logGammaFrac(a,b);
    }

    if( b > 2 )
    {
      if( b < 8 )
      {
        int n = int(b - 1);
        double z = 1.0;
        for( int i=1; i<=n; i++ )
        {
          b -= 1.0;
          z *= b/(a+b);
        }
        return log(z) + log_gamma(a) + log_gamma(b) - logGammaSum(a,b);
      }
      return log_gamma(a) + logGammaFrac(a,b);
    }
    return log_gamma(a) + log_gamma(b) - logGammaSum(a,b);
  }

  if( b >= 8 )
    return log_gamma(a) + logGammaFrac(a,b);

  return log_gamma(a) + log_gamma(b) - log_gamma(a+b);
}

}

}

}

#endif

