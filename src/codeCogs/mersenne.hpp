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
//! Random number generator class using the Mersenne Twister algorithm

#ifndef STATS_RANDOM_MERSENNE_H
#define STATS_RANDOM_MERSENNE_H

#include <assert.h>

#define N           624
#define M           397
#define MATRIX_A    0x9908b0dfUL
#define UPPER_MASK  0x80000000UL
#define LOWER_MASK  0x7fffffffUL
#define MERSENNEDIV 4294967296.0


//! Random number generator class using the Mersenne Twister algorithm

class Mersenne {
  public:

    //! Constructor that initializes the generator with the given seed.

    Mersenne(size_t s = 46875) {
        this->setSeed(s);
        Init((unsigned long)(s));
    }

    //! Copy constructor

    Mersenne(const Mersenne &C) {
        size_t s = C.getSeed();
        this->setSeed (s);
        Init((unsigned long)(s));
    }

    //! Generates an uniform floating point number in the (0, 1) interval (endpoints are excluded).

    double genReal() {
        return ((double)Next() + 0.5) / MERSENNEDIV;
    }

    //! Generates an uniform large integer in the (0, 2^32 - 1) interval.

    unsigned int genInt() {
        return Next();
    }

    //! Resets the seed of the generator to the given parameter.

    void setSeed(size_t s) {
        this->seed_ = s;
        Init((unsigned long)(s));
    }

    //! Returns the seed of the generator.

    size_t getSeed() const {
        return seed_;
    }

  private:
    void Init(unsigned long s){
        mt[0]= s & 0xffffffffUL;
        for (mti = 1; mti < N; ++mti) {
            mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
            mt[mti] &= 0xffffffffUL;
        }
    }

    unsigned long Next(){
        unsigned long y;
        if (mti >= N) {
            int k;
            unsigned long tmp[2] = {0x0UL, MATRIX_A};
            for (k = 0; k < N - M; ++k) {
                y = (mt[k] & UPPER_MASK) | (mt[k + 1] & LOWER_MASK);
                mt[k] = mt[k + M] ^ (y >> 1) ^ tmp[y & 0x1UL];
            }

            for (; k < N - 1; ++k) {
                y = (mt[k] & UPPER_MASK) | (mt[k+1] & LOWER_MASK);
                mt[k] = mt[k + (M - N)] ^ (y >> 1) ^ tmp[y & 0x1UL];
            }

            y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ tmp[y & 0x1UL];
            mti = 0;
        }

        y = mt[mti++];

        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);

        return y;

    }
    size_t seed_;
    unsigned long mt[N];
    int mti;
};


#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK

#endif

