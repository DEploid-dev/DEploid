#include <sstream>      // std::stringstream
#include <stdlib.h>     /* strtol, strtod */
#include <fstream>
#include <stdexcept> // std::invalid_argument
#include <vector>
#include <iostream> // std::cout

#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

#ifndef PARAM
#define PARAM

using namespace std;

class Input{
 friend class McmcSample;
 friend class TestInput;
  public:
    Input( const char plafFileName[],
           const char refFileName[],
           const char altFileName[],
           size_t kStrain );
    ~Input ();

  private:
    // Read in input
    vector <double> plaf;
    vector <double> refCount;
    vector <double> altCount;
    size_t kStrain_;

    void readFileLines(const char inchar[], vector <double> & out_vec);
};


//template <typename T>
//std::ostream& operator<<(std::ostream& os, const std::vector<T> &vec) {
  //typename std::vector<T>::const_iterator it;
  //for (it = vec.begin(); it != vec.end(); ++it) os << *it << " ";
  //return os;
//}

#endif
