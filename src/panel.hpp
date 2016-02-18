#include <vector>

#ifndef PANEL
#define PANEL

using namespace std;

class Panel{
 friend class TestPanel;
 friend class UpdateSingleHap;
 friend class UpdateHap;
  private:
    // Member
    // content is a matrix of n.loci by n.strains, i.e. content length is n.loci
    vector < vector < double > > content_;
    vector < double > recombProbs_;
    size_t nLoci_;
    size_t nPanel_;

  public:
    Panel(const char inchar[]);
    ~Panel(){};

    // Method
    void print();
};

#endif
