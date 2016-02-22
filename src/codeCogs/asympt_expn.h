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
//! Evaluate Del(a) + Del(b) - Del(a+b).


#ifndef MATHS_ALGEBRA_SERIES_ASYMPT_EXPN_H
#define MATHS_ALGEBRA_SERIES_ASYMPT_EXPN_H

#include <assert.h>

#include <math.h>
#include <xsub_ln_add1.h>
#include <errorfnc.h>
#include <errorfnc_exp.h>

namespace Maths
{

namespace Algebra
{

namespace Series
{

//! Evaluate Del(a) + Del(b) - Del(a+b).

double bcorr( double a, double b )
{
  assert(a>=8);
  assert(b>=8);

  const double c0 = 0.0833333333333333;
  const double c1 = -0.00277777777760991;
  const double c2 = 0.000793650666825390;
  const double c3 = -0.000595202931351870;
  const double c4 = 0.000837308034031215;
  const double c5 = -0.00165322962780713;

  double a0,b0;
  if( a < b )
  {
    a0 = a;
    b0 = b;
  }
  else
  {
    a0 = b;
    b0 = a;
  }

  double h = a0/b0;
  double c = h / (1.0+h);
  double x = 1.0 / (1.0+h);
  double x2 = x*x;

  // SET SN = (1 - X**N)/(1 - X)
  double s3 = 1.0 + x + x2;
  double s5 = 1.0 + x + x2*s3;
  double s7 = 1.0 + x + x2*s5;
  double s9 = 1.0 + x + x2*s7;
  double s11 = 1.0 +x + x2*s9;

  // SET W = DEL(B) - DEL(A + B)
  double t = pow( 1.0/b0, 2.0 );
  double w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3)*t + c0;
  w *= (c/b0);

  // DEL(A) + W
  t = pow(1.0e0/a0,2.0);
  return (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/a0 + w;
}

//! Asymptotic Series Expansion for ix(a,b) for large a and b.

double asympt_expn( double a, double b, double lambda, double eps )
{
  const double e0 = 1.12837916709551;
  const double e1 = .353553390593274;
  const int num = 20;
  int K3 = 1;
  double bsum,dsum,f,h,h2,hn,j0,j1,r,r0,r1,s,sum,t,t0,t1,w,w0,z,z0,z2,zn,znm1;
  int i,im1,imj,j,m,mm1,mmj,n,np1;
  double a0[21],b0[21],c[21],d[21],T1,T2;
  double basym = 0.0;
  if( a < b )
  {
    h = a / b;
    r0 = 1.0 / (1.0+h);
    r1 = (b-a) / b;
    w0 = 1.0 / sqrt( a*(1.0+h) );
  }
  else
  {
    h = b / a;
    r0 = 1.0 / (1.0+h);
    r1 = (b-a) / a;
    w0 = 1.0 / sqrt( b*(1.0+h) );
  }
  T1 = -lambda / a;
  T2 = lambda / b;
  f = ( a * Maths::Arithmetic::xsub_ln_add1(T1) )
    + ( b * Maths::Arithmetic::xsub_ln_add1(T2) );
  t = exp(-f);
  if( t == 0)
    return basym;
  else
  {
    z0 = sqrt(f);
    z = 0.5 * z0 / e1;
    z2 = f+f;
    a0[0] = 2.0 / 3.0 * r1;
    c[0] = -0.5 * a0[0];
    d[0] = -c[0];

    if( K3==0 )
      j0 = 0.5 / e0 * Maths::Special::errorFnC( z0 );
    else
      j0 = 0.5 / e0 * Maths::Special::errorFnC_exp( z0 );

    j1 = e1;
    sum = j0 + d[0] * w0 * j1;
    s = 1.0;
    h2 = h*h;
    hn = 1.0e0;
    w = w0;
    znm1 = z;
    zn = z2;
    for( n=2; n<=num; n+=2 )
    {
      hn = h2 * hn;
      a0[n-1] = 2.0 * r0 * (1.0 + h * hn) / ((double)n + 2.0);
      np1 = n+1;
      s += hn;
      a0[np1-1] = 2.0 * r1 * s / ((double)n + 3.0);
      for( i=n; i<=np1; i++ )
      {
        r = -0.5 * ((double)i + 1.0);
        b0[0] = r * a0[0];
        for( m=2; m<=i; m++ )
        {
          bsum = 0.0;
          mm1 = m-1;
          for( j=1; j<=mm1; j++ )
          {
            mmj = m-j;
            bsum += ( (double)j * r - (double)mmj ) * a0[j-1] * b0[mmj-1];
          }
          b0[m-1] = r * a0[m-1] + bsum / (double)m;
        }
        c[i-1] = b0[i-1] / ((double)i + 1.0);
        dsum = 0.0;
        im1 = i-1;
        for( j=1; j<=im1; j++ )
        {
          imj = i-j;
          dsum += d[imj-1] * c[j-1];
        }
        d[i-1] = -dsum - c[i-1];
      }
      j0 = e1 * znm1 + ((double)n - 1.0) * j0;
      j1 = e1 * zn + (double)n * j1;
      znm1 = z2 * znm1;
      zn = z2 * zn;
      w = w0 * w;
      t0 = d[n-1] * w * j0;
      w = w0 * w;
      t1 = d[np1-1] * w * j1;
      sum += t0 + t1;
      if( fabs(t0)+fabs(t1) <= eps*sum )
        break;
    } // end for ( n=2; n<=num; n+=2 )

    return e0 * t * exp( -bcorr(a,b) ) * sum;
  } // end if-else( t == 0)
}

}

}

}

#endif

