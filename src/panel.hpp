#include <vector>
#include <fstream>
#include <stdexcept>      // std::invalid_argument

#ifndef PANEL
#define PANEL

using namespace std;

class Panel{
  private:
    // Member
    // content is a matrix of n.loci by n.strains, i.e. content length is n.loci
    vector < vector < double > > content;

  public:
    Panel(const char inchar[]);
    ~Panel(){};

    // Method
    void print();
};

#endif
