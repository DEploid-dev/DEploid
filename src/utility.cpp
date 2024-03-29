/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of dEploid.
 *
 * dEploid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <iterator>   // std::distance
#include <algorithm>  // find

#include "utility.hpp"
#include "codeCogs/loggammasum.h"  // which includes log_gamma.h
#include "codeCogs/gamma.h"
#include "codeCogs/logbeta.h"


// See http://stackoverflow.com/questions/10847007/
// using-the-gaussian-probability-density-function-in-c
double normal_pdf(double x, double m, double s) {
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-(0.5) * a * a);
}


vector <double> computeCdf(const vector <double> & dist) {
    vector <double> cdf;
    double cumsum = 0;
    for (double p : dist) {
        cumsum += p;
        cdf.push_back(cumsum);
    }
    assert(cdf.size() == dist.size());
    // assert(cdf.back() == 1.0);
    return cdf;
}


double sumOfMat(const vector <vector <double> > & matrix) {
    double tmp = 0.0;
    for (auto const& array : matrix) {
        for (auto const& value : array) {
            tmp += value;
        }
    }
    return tmp;
}


void normalizeBySum(vector <double> & array ) {
    double sumOfArray = sumOfVec(array);
    for (vector<double>::iterator it = array.begin(); it != array.end(); ++it) {
        *it /= sumOfArray;
    }
}


void normalizeByMax(vector <double> & array ) {
    double maxOfArray = max_value(array);
    for (vector<double>::iterator it = array.begin(); it != array.end(); ++it) {
        *it /= maxOfArray;
    }
}


void normalizeBySumMat(vector <vector <double> > & matrix) {
    double tmpsum = sumOfMat(matrix);
    for (size_t i = 0; i < matrix.size(); i++) {
        for (vector<double>::iterator it = matrix[i].begin();
            it != matrix[i].end(); ++it) {
            *it /= tmpsum;
        }
    }
}


vector <double> calcLLKs(const vector <double> &refCount,
                         const vector <double> &altCount,
                         const vector <double> &expectedWsaf,
                         size_t firstIndex, size_t length,
                         double fac, double err) {
    assert(length <= expectedWsaf.size());
    vector <double> tmpLLKs(length, 0.0);
    size_t index = firstIndex;
    for (size_t i = 0; i < length; i++) {
        assert(expectedWsaf[i] >= 0);
        // assert (expectedWsaf[i] <= 1);
        tmpLLKs[i] = log(calcSiteLikelihood(refCount[index], altCount[index],
                                            expectedWsaf[i], err, fac));
        index++;
    }
    return tmpLLKs;
}

vector <log_double_t> calcSiteLikelihoods(const vector <double> &refCount,
                                          const vector <double> &altCount,
                                          const vector <double> &expectedWsaf,
                                          size_t firstIndex, size_t length,
                                          double fac, double err) {
    assert(expectedWsaf.size() == length);
    vector <log_double_t> siteLikelihoods(length);
    size_t index = firstIndex;
    for (size_t i = 0; i < length; i++) {
        assert(expectedWsaf[i] >= 0);
        // assert (expectedWsaf[i] <= 1);
        siteLikelihoods[i] = calcSiteLikelihood(refCount[index], altCount[index],
                                                expectedWsaf[i], err, fac);
        index++;
    }
    return siteLikelihoods;
}

log_double_t Beta(double x, double y)
{
    return exp_to<log_double_t>(Maths::Special::Gamma::logBeta(x,y));
}

// Probability of observing k from a beta-binominal(n,a,b)
log_double_t beta_binomial_pr(int n, double a, double b, int k)
{
    if (k < 0) return 0;
    if (k > n) return 0;

    // pr = choose(n,k) * beta(k+a, n-k+b) / beta(a,b);
    auto pr = Beta(k+a, n-k+b) / Beta(a,b);

    // choose(n,k) = 1/[(n+1) * beta(n-k+1,k+1)]
    pr /= Beta(n-k+1, k+1);
    pr /= (n+1);

    // Previously we were ignoring the choose(n,k) term.
    // This might have been OK, since probabilities of (a,b) would be correct up to a constant.
    // But lets include it, to be safe.

    return pr;
}

log_double_t calcSiteLikelihood(double ref, double alt, double unadjustedWsaf, double err,
    double fac) {

    // Adjusting for sequencing error
    double adjustedWsaf = unadjustedWsaf+err*(1-2*unadjustedWsaf);

    return beta_binomial_pr(ref + alt, adjustedWsaf*fac, (1-adjustedWsaf)*fac, alt);
}


size_t sampleIndexGivenProp(RandomGenerator* rg, vector <double> proportion) {
    #ifndef NDEBUG
        auto biggest = std::max_element(
            std::begin(proportion), std::end(proportion));
        return std::distance(proportion.begin(), biggest);
    #else
        vector <size_t> tmpIndex;
        for ( size_t i = 0; i < proportion.size(); i++ ) {
            tmpIndex.push_back(i);
        }
        vector <double> tmpCdf = computeCdf(proportion);

        double u = rg->sample();
        size_t i = 0;
        for ( ; i < tmpCdf.size() ; i++) {
            if ( u < tmpCdf[i] ) {
                break;
            }
        }
        return i;
    #endif
}


vector <double> reshapeMatToVec(const vector < vector <double> > &Mat) {
    vector <double> tmp;
    for (auto const& array : Mat) {
        for (auto const& value : array) {
            tmp.push_back(value);
        }
    }
    return tmp;
}


double betaPdf(double x, double a, double b) {
    assert(x >= 0 && x <= 1);
    assert(a >= 0);
    assert(b >= 0);
    double p = Maths::Special::Gamma::gamma(a + b) /
        (Maths::Special::Gamma::gamma(a) * Maths::Special::Gamma::gamma(b));
    double q = pow(1 - x, b - 1) * pow(x, a - 1);
    return p * q;
}


double logBetaPdf(double x, double a, double b) {
    assert(x >= 0 && x <= 1);
    assert(a >= 0);
    assert(b >= 0);
    double ret = Maths::Special::Gamma::log_gamma(a+b) -
                 Maths::Special::Gamma::log_gamma(a) -
                 Maths::Special::Gamma::log_gamma(b) +
                 (b-1) * log(1-x) + (a-1) * log(x);
    return ret;
}


double binomialPdf(int s, int n, double p) {
    assert(p >= 0 && p <= 1);
    double ret = n_choose_k(n, s);
    ret *= pow(p, static_cast<double>(s));
    ret *= pow((1.0-p), static_cast<double>(n-s));
    return ret;
}


// double betaDistConst( double a , double b) {
    // double ret = Maths::Special::Gamma::gamma(a + b) /
    // (Maths::Special::Gamma::gamma(a) * Maths::Special::Gamma::gamma(b));
    // return ret;
// }


double rBeta(double alpha, double beta, RandomGenerator* rg) {
    double mxAt = (alpha-1.0)/(alpha+beta-2.0);
    double mx = betaPdf(mxAt, alpha, beta);
    double y, u;
    do {
        u = rg->sample();
        y = rg->sample();
    } while ( u > (betaPdf(y, alpha, beta)/mx));

    return y;
}
