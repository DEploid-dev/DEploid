#include <vector>
#include <iostream>
#include "utility.hpp"
#include "panel.hpp"
#include "mersenne_twister.hpp"

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
               vector < vector <double> > &haplotypes, MersenneTwister* rg, Panel* panel);
    virtual ~UpdateHap(){}

    Panel* panel_;
    double missCopyProb;
    MersenneTwister* rg_;
    size_t strainIndex_;
    size_t kStrain_;
    size_t nLoci_;
    size_t nPanel_;

    vector < vector <double> > emission_;
    vector < vector <double> > fwdProbs_;

    // Methods
    virtual void findUpdatingStrain( vector <double> proportion ){};
    virtual void calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes){};
    virtual void calcHapLLKs( vector <double> &refCount, vector <double> &altCount){};
    virtual void buildEmission(){};
    virtual void calcFwdProbs(){};
    virtual void samplePaths(){};
    virtual void addMissCopying(){};
};


class UpdateSingleHap : public UpdateHap{
  public:
     UpdateSingleHap( vector <double> &refCount,
                      vector <double> &altCount,
                      vector <double> &expectedWsaf,
                      vector <double> &proportion,
                      vector < vector <double> > &haplotypes, MersenneTwister* rg, Panel* panel);
    ~UpdateSingleHap(){}
  private:
    size_t strainIndex_;
    vector <double> expectedWsaf0_;
    vector <double> expectedWsaf1_;
    vector <double> llk0_;
    vector <double> llk1_;

    vector <double> path_;
    vector <double> hap_;
    vector <double> newLLK;

    // Methods
    void findUpdatingStrain( vector <double> proportion );
    void calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes);
    void calcHapLLKs( vector <double> &refCount, vector <double> &altCount);
    void buildEmission();
    void calcFwdProbs();
    void samplePaths();
    void addMissCopying();
};



#endif
