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
//! Evaluates a polynomial of degree N.


#ifndef MATHS_ALGEBRA_POLYNOMIAL_POLY_EVAL_H
#define MATHS_ALGEBRA_POLYNOMIAL_POLY_EVAL_H

namespace Maths
{

namespace Algebra
{

namespace Polynomial
{

//! Evaluates a polynomial of degree N.

double polyEval(double x, const double coef[], int N)
{
  const double *p = coef;
  double ans = *p++;
  int i = N;

  do
    ans = ans * x  +  *p++;
  while( --i );

  return( ans );
}

//! Evaluates a polynomial of degree N, with $C_N=1.0$.

double polyEval1(double x, const double coef[], int N)
{
  const double *p = coef;
  double ans = x + *p++;
  int i = N-1;

  do
    ans = ans * x  + *p++;
  while( --i );

  return( ans );
}

}

}

}

#endif

