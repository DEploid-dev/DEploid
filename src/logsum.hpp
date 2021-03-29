#ifndef LOGSUM_H
#define LOGSUM_H

#include <cmath>
#include <cassert>

constexpr double max_float = 3.40282347e+38F;
constexpr double log_limit = -max_float/100;
constexpr double log_0 = -max_float;

// NATS is 52*log(2) for 52 bits of precision
// HMM... the long doubles have 64 bits of precision...
constexpr double NATS = 40;

inline double logsum_nocheck(double x, double y) {
    if (std::abs(x-y) > NATS)
	return ((x > y) ? x : y);
    else
	return (x + log1p(exp(y - x)));
}

inline double logsum(double x, double y)
{
    double temp = y-x;
    if (temp > NATS or x < log_limit)
	return y;
    else if (temp < -NATS or y < log_limit)
	return x;
    else
	return (x + log1p(exp(temp)));
}

inline void loginc(double& x, double y)
{
    double temp = y-x;
    if (temp > NATS or x < log_limit)
	x=y;
    else if (temp < -NATS or y < log_limit)
	;
    else
	x += log1p(exp(temp));
}

inline double logdiff(double x, double y) {
    assert(x >= y);
    double temp = y-x;
    if (temp < -NATS or y < log_limit)
	return x;
    else if (x == y)
	return log_0;
    else
	return (x + log1p(-exp(temp)));
}
#endif
