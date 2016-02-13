#include "param.hpp"

Input::Input( const char plafFileName[],
              const char refFileName[],
              const char altFileName[],
              size_t kStrain ){
    (void)readFileLines( refFileName, this->refCount);
    //dout <<" refCount.size() = "<< refCount.size()<<endl;
    (void)readFileLines( altFileName, this->altCount);
    //dout <<" altCount.size() = "<< altCount.size()<<endl;
    (void)readFileLines( plafFileName, this->plaf);
    //dout <<" plaf.size() = "<< plaf.size()<<endl;
    this->kStrain_ = kStrain;
}


Input::~Input(){}


void Input::readFileLines(const char inchar[], vector <double> & out_vec){
    ifstream in_file( inchar );
    string tmp_line;
    if ( in_file.good() ){
        getline ( in_file, tmp_line ); // skip the first line, which is the header
        getline ( in_file, tmp_line );
        while ( tmp_line.size() > 0 ){
            size_t field_start = 0;
            size_t field_end = 0;
            size_t field_index = 0;

            while ( field_end < tmp_line.size() ){
                field_end = min ( tmp_line.find('\t',field_start), tmp_line.find('\n', field_start) );
                string tmp_str = tmp_line.substr( field_start, field_end - field_start );

                if ( field_index == 2 ){
                    out_vec.push_back( strtod(tmp_str.c_str(), NULL) );
                }
                field_start = field_end+1;
                field_index++;
          }
          getline ( in_file, tmp_line );
      }
    } else {
        throw std::invalid_argument("Invalid input file. " + string (inchar) );

    }
    in_file.close();
}
