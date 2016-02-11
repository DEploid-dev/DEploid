#include <vector>
#include <iostream>
#include "mersenne_twister.hpp"


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
  public:
    McmcSample( size_t nSample = 100, size_t mcmcSampleRate = 5, size_t kStrain = 5, size_t seed = 88 );
    ~McmcSample();
    void runMcmcChain( );

  private:
    // Variables
    size_t kStrain_;

    double burnIn_;
    size_t maxIteration_;
    size_t mcmcThresh_;
    size_t mcmcSampleRate_;
    int eventInt_;
    MersenneTwister* rg_;

    //EventType currentUpdateEvent;
    size_t seed_;

    size_t currentMcmcIteration_;

    vector <double> currentTitre_;
    double currentLodPriorTitre_;
    vector <double> currentProp_;
    vector < vector <double> > currentHap;

    // Methods
    void sampleMcmcEvent();
    void recordMcmcSample();
    void initializeProp();
    void initializeTitre();

    void calcMaxIteration( size_t nSample, size_t mcmcSampleRate );

    // Moves
    void updateProportion();
    void updateSingleHap();
    void updatePairHaps();


//fun.initialiseBundle <- function ( Data, PLAF, n.event, initialK = 3, labelP = LABELP){
  //n.loci<-length(PLAF);

  //#Initialise titres
  //titre<-rnorm(initial.k, MN_LOG_TITRE, SD_LOG_TITRE);
//###  titre<-rep(0, initial.k) # Set equal titres
  //initial.propotion<-fun.titre2prop(titre);
  //log.prior.titre<-sum(dnorm(titre, MN_LOG_TITRE, SD_LOG_TITRE, log=TRUE));


//#  initial.haplotype<-array(0, c(n.loci, K_MAX));
  //initial.haplotype<-array(0, c(n.loci, initial.k));
  //initial.haplotype_counter<-array(0, c((n.loci+1), initial.k));
  //initial.path_counter<-array(0, c((n.loci+1), initial.k));

  //if (initial.k>0) for (i in 1:initial.k) initial.haplotype[,i]<-rbinom(n.loci, 1, PLAF);

  //initial.expected.WSAF<-fun.calc.f.samp(initial.haplotype, initial.propotion);
//#  initial.likelihood<-fun.llk.with.error(Data, initial.expected.WSAF, Beta.LLK.with.u,  type=mode);
  //initial.likelihood<-fun.llk( Data, initial.expected.WSAF, type = mode)

  //initial.accpt.rate<-rep(0, n.event);


  //bundle = list( k = initial.k,                        #1
                 //titre = titre,                        #2
                 //prop = initial.propotion,             #3
                 //h = initial.haplotype,                #4
                 //llk.old = initial.likelihood,         #5
                 //log.prior.titre = log.prior.titre,    #6
                 //accpt.rate = initial.accpt.rate,      #7
                 //expected.WSAF = initial.expected.WSAF,#8
                 //count = initial.haplotype_counter,    #9
                 //path = initial.path_counter,          #10
                 //switch.updateOne.one = NA,
                 //switch.updateTwo.one = NA,
                 //switch.updateTwo.two = NA
               //)
  //return(bundle)
//}

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

#endif

