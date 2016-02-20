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
//! Returns gamma function of the argument.

#ifndef MATHS_SPECIAL_GAMMA_GAMMA_H
#define MATHS_SPECIAL_GAMMA_GAMMA_H

#ifndef cc_error
#include <assert.h>
#ifndef NDEBUG
#include <iostream>
#define cc_error(a) std::cerr << (a) << std::endl;
#else
#define cc_error(a)
#endif
#endif

#include <poly_eval.h>
#include <stirling.h>
#include <math.h>
#include <float.h>

#define MAXNUM 1.79769313486231570815E308
#define PI 3.141592653589793238

namespace Maths
{

namespace Special
{

namespace Gamma
{

//! Returns gamma function of the argument.

double gamma(double x, int *sign=NULL)
{
  static double P[] = {
    1.60119522476751861407E-4,
    1.19135147006586384913E-3,
    1.04213797561761569935E-2,
    4.76367800457137231464E-2,
    2.07448227648435975150E-1,
    4.94214826801497100753E-1,
    9.99999999999999996796E-1
  };
  static double Q[] = {
    -2.31581873324120129819E-5,
    5.39605580493303397842E-4,
    -4.45641913851797240494E-3,
    1.18139785222060435552E-2,
    3.58236398605498653373E-2,
    -2.34591795718243348568E-1,
    7.14304917030273074085E-2,
    1.00000000000000000320E0
  };

  double p, q, z;
  int i;

  int isign= 1;
  if(sign!=NULL) *sign=isign;
#ifdef NANS
    if( isnan(x) )
      return(x);
#endif
#ifdef INFINITIES
#ifdef NANS
    if( x == INFINITY )
      return(x);
    if( x == -INFINITY )
      return(NAN);
#else
    if( !isfinite(x) )
      return(x);
#endif
#endif
    q = fabs(x);

    if( q > 33.0 )
    {
      if( x < 0.0 )
      {
        p = floor(q);
        if( p == q )
        {
#ifdef NANS
gamnan:
          cc_error("gamma(): Nan ");
          return (NAN);
#else
          goto goverf;
#endif
        }
        i = (int)p;
        if( (i & 1) == 0 )
          isign = -1;

        if(sign!=NULL)
          *sign=isign;

        z = q - p;
        if( z > 0.5 )
        {
          p += 1.0;
          z = q - p;
        }
        z = q * sin( PI * z );
        if( z == 0.0 )
        {
#ifdef INFINITIES
          return( isign * INFINITY);
#else
goverf:
          cc_error("gamma(): overflow ");

          return( isign * MAXNUM);
#endif
        }
        z = fabs(z);
        z = PI/(z * Maths::Special::Gamma::stirling(q) );
      }
      else
      {
        z = Maths::Special::Gamma::stirling(x);
      }
      return( isign * z );
    }

    z = 1.0;
    while( x >= 3.0 )
    {
      x -= 1.0;
      z *= x;
    }

    while( x < 0.0 )
    {
      if( x > -1.E-9 )
        goto jmpsmall;
      z /= x;
      x += 1.0;
    }

    while( x < 2.0 )
    {
      if( x < 1.e-9 )
        goto jmpsmall;
      z /= x;
      x += 1.0;
    }

    if( x == 2.0 )
      return(z);

    x -= 2.0;
    p = Maths::Algebra::Polynomial::polyEval( x, P, 6 );
    q = Maths::Algebra::Polynomial::polyEval( x, Q, 7 );
    return( z * p / q );

jmpsmall:
      if( x == 0.0 )
      {
#ifdef INFINITIES
#ifdef NANS
        goto gamnan;
#else
        return( INFINITY );
#endif
#else
        cc_error("gamma(): singularity ");

        return( MAXNUM );
#endif
      }
    else
      return( z/((1.0 + 0.5772156649015329 * x) * x) );
}

}

}

}

#undef MAXNUM
#undef PI

#endif

