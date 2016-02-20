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

#ifndef MATHS_SPECIAL_ERRORFN_H
#define MATHS_SPECIAL_ERRORFN_H

#include <assert.h>
#include <poly_eval.h>

namespace Maths
{

namespace Special
{

//! The Error Function.

double errorFn( double x )
{
  assert(x>-3 && x<3);
  static double T[] = {
   9.60497373987051638749E0,
   9.00260197203842689217E1,
   2.23200534594684319226E3,
   7.00332514112805075473E3,
   5.55923013010394962768E4
  };
  static double U[] = {
  /* 1.00000000000000000000E0,*/
   3.35617141647503099647E1,
   5.21357949780152679795E2,
   4.59432382970980127987E3,
   2.26290000613890934246E4,
   4.92673942608635921086E4
  };
  double z = x * x,
  y = x * Maths::Algebra::Polynomial::polyEval( z, T, 4 )
        / Maths::Algebra::Polynomial::polyEval1( z, U, 5 );
  return( y );
}

}

}

#endif

