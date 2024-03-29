#ifndef LOG_DOUBLE_H
#define LOG_DOUBLE_H

#include <cassert>
#include <iostream>
#include "logsum.hpp"

/// A class for handling positive real number in terms of their natural log.
class log_double_t {
    /// Natural log of the number.
    double value;
public:

    /// Access the log of the number.
    constexpr double  log() const {return value;}
    /// Access the log of the number.
    constexpr double& log()       {return value;}

    log_double_t& operator +=(const log_double_t& y) {loginc(value,y.log()); return *this;}
    log_double_t& operator -=(const log_double_t& y) {value = logdiff(log(),y.log()); return *this;}
    constexpr log_double_t& operator *=(const log_double_t& y) {value += y.log(); return *this;}
    constexpr log_double_t& operator /=(const log_double_t& y) {value -= y.log(); return *this;}

    operator double() const {return exp(value);}

    constexpr log_double_t& operator=(double y) {
	assert(y >= 0);
	if (y == 0)
	    value = log_0;
	else if (y == 1)
	    value = 0;
	else
	    value = ::log(y);
	return *this;
    }

    constexpr log_double_t():value(log_0) {}

    constexpr log_double_t(double x):value(0) {
	operator=(x);
    }
};

constexpr double log(log_double_t x) {
    return x.log();
}

#define decl_double(rtype,op)					\
    inline rtype operator op (double x,log_double_t y) {	\
	return log_double_t(x) op y;				\
    }								\
								\
    inline rtype operator op (log_double_t x,double y) {	\
	return x op log_double_t(y);				\
    }

inline log_double_t operator+(log_double_t x,log_double_t y) {
    log_double_t z = x;
    z += y;
    return z;
}

inline log_double_t operator-(log_double_t x,log_double_t y) {
    log_double_t z = x;
    z -= y;
    return z;
}

constexpr log_double_t operator*(log_double_t x,log_double_t y) {
    log_double_t z = x;
    z *= y;
    return z;
}

constexpr log_double_t operator/(log_double_t x,log_double_t y) {
    log_double_t z=x;
    z /= y;
    return z;
}

decl_double(log_double_t,+)
decl_double(log_double_t,-)
decl_double(log_double_t,*)
decl_double(log_double_t,/)

constexpr bool operator< (log_double_t x,log_double_t y) {
    return log(x)<log(y);
}

decl_double(bool,<)

constexpr bool operator<=(log_double_t x,log_double_t y) {
    return log(x)<=log(y);
}

decl_double(bool,<=)

constexpr bool operator> (log_double_t x,log_double_t y) {
    return log(x)>log(y);
}

constexpr bool operator>(log_double_t x,double y) {
    if (y==0)
	return log(x) > log_limit;
    else
	return x > log_double_t(y);
}

constexpr bool operator>(double  x,log_double_t y) {
    return log_double_t(x) > y;
}

constexpr bool operator>=(log_double_t x,log_double_t y) {
    return log(x)>=log(y);
}

decl_double(bool,>=)


constexpr bool operator==(log_double_t x,log_double_t y) {
    return log(x)==log(y);
}

inline bool operator==(log_double_t x,double y) {
    if (y==0.0)
	return log(x) <= log_limit;
    else if (y==1)
	return log(x) == 0.0;
    else
	return log(x) == log(y);
}

/*
  constexpr bool operator==(log_double_t x,double y) {
  if (y==0)
  return (x==log_0);
  else if (y==1)
  return (x==0);
  else
  return (x==log(y));
  }
*/

constexpr bool operator!=(log_double_t x,log_double_t y) {
    return log(x)!=log(y);
}

inline bool operator!=(log_double_t x,double y) {
    if (y==0.0)
	return log(x) > log_limit;
    else if (y==1)
	return log(x) != 0.0;
    else
	return log(x) != log(y);
}

#undef decl_double

constexpr log_double_t pow(log_double_t x,double p) 
{
    x.log() *= p;
    return x;
}

inline bool different(log_double_t x,log_double_t y,double tol=1.0e-9)
{
    double diff = log(x) - log(y);

    if (std::abs(diff)<tol) return false;

    if (x == y) return false;

    return true;
}

template<class T> T exp_to(double);

template<> constexpr log_double_t exp_to<log_double_t>(double x) {
    log_double_t y;
    y.log() = x;
    return y;
}

// Don't do pow<log_double_t>(x,y): do pow(log_double_t(x),y), instead.
//
// template<class T> T pow(double,double);
// 
// template<> constexpr log_double_t pow<log_double_t>(double x,double p) {
//   log_double_t y;
//   y.log() = p * ::log(x);
//   return y;
// }

inline std::ostream& operator<<(std::ostream& o,const log_double_t& e) {
    return o<<log(e);
}

#endif
