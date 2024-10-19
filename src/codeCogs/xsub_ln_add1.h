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
//! Evaluation of the function x - ln(1 + x)

#ifndef MATHS_ARITHMETIC_XSUB_LN_ADD1_H
#define MATHS_ARITHMETIC_XSUB_LN_ADD1_H

#include <math.h>

namespace Maths
{

namespace Arithmetic
{

//! Evaluation of the function x - ln(1 + x)

double xsub_ln_add1( double x )
{
  if( x < -0.39 || x > 0.57 )
    return x - log( x + 1.0 );

  double h, w1;

  if( x < -0.18 )
  {
    const double a = 0.0566749439387324;
    h = (x + 0.3) / 0.7;
    w1 = a - h * 0.3;
  }
  else
  {
    if( x > 0.18 )
    {
      const double b = 0.456512608815524;
      h = 0.75 * x - 0.25;
      w1 = b + h / 3.0;
    }
    else
    {
      h = x;
      w1 = 0.0;
    }
  }
  const double p0 = 0.333333333333333;
  const double p1 = -0.224696413112536;
  const double p2 = 0.00620886815375787;
  const double q1 = -1.27408923933623;
  const double q2 = 0.354508718369557;

  double r = h / (h+2);
  double t = r*r;
  double w = (  ( p2*t + p1 )*t + p0  ) / (  ( q2*t + q1 )*t + 1.0  );

  return 2.0*t*( 1.0/(1.0-r) - r*w ) + w1;
}

}

}

#endif

