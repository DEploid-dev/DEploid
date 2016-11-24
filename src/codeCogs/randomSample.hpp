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
//! Generates random numbers following a standard normal distribution.

#ifndef STATS_DISTS_CONTINUOUS_NORMAL_RANDOMSAMPLE_H
#define STATS_DISTS_CONTINUOUS_NORMAL_RANDOMSAMPLE_H

#include <mersenne.hpp>
#include <math.h>

#define NORMALDENS(x) ((fabs(x) > 8.0) ? 0 : 0.398942280 * exp(- x * x / 2))


//! Generates random numbers following a standard normal distribution.

class StandNormalRandomSample : public Mersenne {
  public:

    //! Class constructor.
    StandNormalRandomSample(size_t s = 65234) : Mersenne(s) {
        sx = new double[60];
        sfx = new double[60];
        double sxi = 0.0;
        double inc = 0.01;
        for (int i = 0; i < 60; ++i) {

            sx[i] = sxi;
            double f1 = NORMALDENS(sxi);
            sfx[i] = f1;
            if (f1 <= 0.0) {
                xi = 2 * i;
                return;
            }
            sxi += inc / f1;
        }
    }

    //! Class destructor.
    ~StandNormalRandomSample() {
        delete [] sx;
        delete [] sfx;
    }

    //! Generates a random deviate from the standard normal distribution.
    double genReal() {
        int ir;
        double s, ak, y;
        do {

              s = 1.0;
              double r1 = Mersenne::genReal();
              if (r1 > 0.5) {
                 s = -1.0;
                 r1 = 1.0 - r1;
              }
              ir = (int)(r1 * xi);
              double sxi = sx[ir];
              ak = sxi + (sx[ir + 1] - sxi) * Mersenne::genReal();
              y = sfx[ir] * Mersenne::genReal();

        } while (y >= sfx[ir + 1] && y >= NORMALDENS(ak));

        return s * ak;
    }

  private:
    double xi, *sx, *sfx;
};


#endif

