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
//! The Complementary Error Function.

#ifndef MATHS_SPECIAL_ERRORFNC_H
#define MATHS_SPECIAL_ERRORFNC_H

#ifndef cc_error
#include <assert.h>
#ifndef NDEBUG
#include <iostream>
#define cc_error(a) std::cerr << (a) << std::endl;
#else
#define cc_error(a)
#endif
#endif

#include <expx2.h>
#include <errorfnc_exp.h>
#include <errorfn.h>

#define MINLOG -7.08396418532264106224E2

namespace Maths
{

namespace Special
{

//! The Complementary Error Function.

double errorFnC(double x)
{
  double a;
  if( x < 0 ) a = -x;
  else a = x;
  if( a < 1.0 )
    return( 1.0 - Maths::Special::errorFn(x) );

  double z = -x * x;
  if( z < MINLOG )
  {
under:
    cc_error( "errorFnC: UNDERFLOW" );
    if( x < 0 )
      return( 2 );
    else
      return( 0 );
  }

  z = Maths::Arithmetic::expx2( x, -1 );

  double y = Maths::Special::errorFnC_exp(a) * z;
  if( x < 0 )
    y = 2.0 - y;

  if( y == 0 )
    goto under;

  return(y);
}

}

}

#undef MINLOG

#endif

