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
//! The Exponential Complementary Error Function.

#ifndef MATHS_SPECIAL_ERRORFNC_EXP_H
#define MATHS_SPECIAL_ERRORFNC_EXP_H

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

namespace Maths
{

namespace Special
{

//! The Exponential Complementary Error Function.

inline double errorFnC_exp(double x)
{
  static double P[] = {
   2.46196981473530512524E-10,
   5.64189564831068821977E-1,
   7.46321056442269912687E0,
   4.86371970985681366614E1,
   1.96520832956077098242E2,
   5.26445194995477358631E2,
   9.34528527171957607540E2,
   1.02755188689515710272E3,
   5.57535335369399327526E2
  };
  static double Q[] = {
  /* 1.00000000000000000000E0,*/
   1.32281951154744992508E1,
   8.67072140885989742329E1,
   3.54937778887819891062E2,
   9.75708501743205489753E2,
   1.82390916687909736289E3,
   2.24633760818710981792E3,
   1.65666309194161350182E3,
   5.57535340817727675546E2
  };
  static double R[] = {
   5.64189583547755073984E-1,
   1.27536670759978104416E0,
   5.01905042251180477414E0,
   6.16021097993053585195E0,
   7.40974269950448939160E0,
   2.97886665372100240670E0
  };
  static double S[] = {
  /* 1.00000000000000000000E0,*/
   2.26052863220117276590E0,
   9.39603524938001434673E0,
   1.20489539808096656605E1,
   1.70814450747565897222E1,
   9.60896809063285878198E0,
   3.36907645100081516050E0
  };
  double p,q;
  if( x < 8.0 )
  {
    p = Maths::Algebra::Polynomial::polyEval( x, P, 8 );
    q = Maths::Algebra::Polynomial::polyEval1( x, Q, 8 );
  }
  else
  {
    p = Maths::Algebra::Polynomial::polyEval( x, R, 5 );
    q = Maths::Algebra::Polynomial::polyEval1( x, S, 6 );
  }
  return (p/q);
}

}

}

#endif

