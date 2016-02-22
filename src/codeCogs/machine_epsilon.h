// CodeCogs Commercial License Agreement
// Copyright (C) 2004-2010 CodeCogs, Zyba Ltd, Broadwood, Holford, TA5 1DU, England.
//
// This software is licensed to Joe Zhu 
// for commercial usage by version 1.2.1 of the CodeCogs Commercial Licence. You must 
// read this License (available at www.codecogs.com) before using this software.
//
// If you distribute this file it is YOUR responsibility to ensure that all 
// recipients have a valid number of commercial licenses. You must retain a
// copy of this licence in all copies you make. 
//
// This program is distributed WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the CodeCogs Commercial Licence for more details.
//---------------------------------------------------------------------------------

#ifndef COMPUTING_LOWLEVEL_MACHINE_EPSILON_H
#define COMPUTING_LOWLEVEL_MACHINE_EPSILON_H

namespace Computing
{

namespace Lowlevel
{

//! Returns the round off unit for double precision arithmetic.

inline double machine_epsilon()
{
  double eps = 1;
  while (1 + (eps /= 2) > 1);
  return 2 * eps;
}

}

}

#endif

