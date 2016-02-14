#include <vector>
#include <iostream>
#include "mersenne_twister.hpp"
#include "param.hpp"

#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

#ifndef MCMC
#define MCMC

using namespace std;

//enum EventType {PROPORTION, SINGLE, PAIR};


class McmcSample {
  friend class McmcMachinery;
  private:
    vector< vector <double> > proportion;
    vector< double > sumLLKs;
};

class McmcMachinery {
  public:
    McmcMachinery( Input* input,
                size_t nSample = 100, size_t McmcMachineryRate = 5, size_t seed = 88 );
    ~McmcMachinery();
    void runMcmcChain( );

  private:
  /* Variables */
    Input* input_;
    size_t kStrain_;
    size_t nLoci_;

    double burnIn_;
    size_t maxIteration_;
    size_t mcmcThresh_;
    size_t McmcMachineryRate_;
    int eventInt_;

    size_t seed_;
    MersenneTwister* rg_;

    //EventType currentUpdateEvent;

    size_t currentMcmcIteration_;

    vector <double> currentTitre_;
    double currentLogPriorTitre_;
    vector <double> currentProp_;
    vector <double> currentLLks_;
    vector < vector <double> > currentHap_;
    vector < double > currentExpectedWsaf_;

    std::default_random_engine* std_generator_;// (this->seed_);

    std::normal_distribution<double>* std_normal_distribution_;// (MN_LOG_TITRE, SD_LOG_TITRE);

  /* Methods */
   /* Initialize */
    void initializeMcmcChain();
     void initializeProp();
     void initializeTitre();
     void initializeHap();
     void initializellk();
     void initializeExpectedWsaf();

    void calcExpectedWsaf();
    void calcCurrentLLKs();
    double rBernoulli(double p);
    double calcLLK( double ref, double alt, double unadjustedWsaf, double err = 0.01, double fac=100 ) ;

    void sampleMcmcEvent();
    void recordMcmcMachinery();

    void calcMaxIteration( size_t nSample, size_t McmcMachineryRate );

    double sum( vector <double>& array ){
        double tmp = 0.0;
        for (auto const& value: array){
            tmp += value;
        }
        return tmp;
    }

    void normalizeBySum ( vector <double> & array ){
        double sumOfArray = sum(array);
        for( vector<double>::iterator it = array.begin(); it != array.end(); ++it) {
            *it /= sumOfArray;
        }
    }

    void printArray ( vector <double> array ){
        for (auto const& value: array){
            cout << value << " ";
        }
        cout << endl;
    }

    double MN_LOG_TITRE;
    double SD_LOG_TITRE;

    void calLogPriorTitre();

  /* Moves */
    void updateProportion();
    void updateSingleHap();
    void updatePairHaps();

};


//class BasicData{
    //vector <double> chrom;
    //vector <double> pos;
    //vector <double> value;
//};


//class Falciparum{
    //double prop;
    //vector <unsigned short> hap;
//};

//class MixedFalciparum{
    //vector <Falciparum> falciparums;

//};

//class McmcEvent{

//};

template <typename T> // See http://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
T normal_pdf(T x, T m, T s)
{
    static const T inv_sqrt_2pi = 0.3989422804014327;
    T a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-T(0.5) * a * a);
}

#endif

