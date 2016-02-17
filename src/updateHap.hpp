#include <vector>
#include <iostream>
#include "utility.hpp"
#include "mersenne_twister.hpp"
//#include "param.hpp"

#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

#ifndef HAP
#define HAP

using namespace std;


class UpdateHap{
  friend class McmcMachinery;
  friend class UpdateSingleHap;
    UpdateHap( vector <double> &refCount,
               vector <double> &altCount,
               vector <double> &expectedWsaf,
               vector <double> &proportion,
               vector < vector <double> > &haplotypes, MersenneTwister* rg);
    virtual ~UpdateHap(){}

  //private:
    double missCopyProb;
    MersenneTwister* rg_;
    size_t strainIndex_;
    size_t kStrain_;
    size_t nLoci_;

    vector < vector <double> > emission_;

    // Methods
    virtual void buildEmission(){};
    virtual void calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes){};
    virtual void findUpdatingStrain(){};

virtual void calcHapLLKs( vector <double> &refCount,
                                vector <double> &altCount){};
};


class UpdateSingleHap : public UpdateHap{
  public:
     UpdateSingleHap( vector <double> &refCount,
                      vector <double> &altCount,
                      vector <double> &expectedWsaf,
                      vector <double> &proportion,
                      vector < vector <double> > &haplotypes, MersenneTwister* rg);
    ~UpdateSingleHap(){}
  private:
    vector <double> expectedWsaf0_;
    vector <double> expectedWsaf1_;
    vector <double> llk0_;
    vector <double> llk1_;

    void buildEmission();
    void calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes);
    void findUpdatingStrain();

    void calcHapLLKs( vector <double> &refCount, vector <double> &altCount);
};



#endif
