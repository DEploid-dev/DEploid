
#include <vector>

using namespace std;

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

class McmcMachinary{

};

//enum EventType {PROPORTION, SINGLE, PAIR};

class McmcSample {
  public:
    McmcSample( size_t kStrain = 5, size_t seed = 88 ){
        this->kStrain = kStrain;
        this->seed = seed;
        rg = new MersenneTwister(this->seed);
    };

    ~McmcSample(){
        delete rg;
    }

    void sampleMcmcEvent(){
        this->eventInt = this->rg->sampleInt(3);
        if ( this->eventInt == 0 ){
            this->updateProportion();
        } else if ( this->eventInt == 1 ){
            this->updateSingleHap();
        } else if ( this->eventInt == 2 ){
            this->updatePairHaps();
        } else {
            throw(" should never reach here!!!");
            cout << "eventInt = " << this->eventInt << endl;
        }
    }


    void recordMcmcSample(){
        cout << "Record mcmc sample " <<endl;
    }


  private:
    // Variables
    int eventInt;
    MersenneTwister* rg;
    //EventType currentUpdateEvent;
    size_t seed;
    size_t kStrain;
    vector <double> currentTitre;
    double currentLodPriorTitre;
    vector <double> currentProp;
    vector < vector <double> > currentHap;

    // Methods
    void initializeProp( size_t kStrain ){
        //initial.propotion[1:initial.k]<-fun.sim_prop(initial.k, TRUE);
        //initial.propotion<-initial.propotion/sum(initial.propotion);
    }

    void initializeTitre( size_t kStrain ){

    }

    void updateProportion(){
        cout << "update Proportion "<<endl;
    }

    void updateSingleHap(){
        cout << "update Single Hap "<<endl;
    }

    void updatePairHaps(){
        cout << "update Pair Hap "<<endl;
    }
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
