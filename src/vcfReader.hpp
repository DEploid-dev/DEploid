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

#include <stdlib.h>     /* strtol, strtod */
#include <string>  /* string */
#include <vector>  /* vector */
#include <fstream>
#include "exceptions.hpp"
#include "variantIndex.hpp"
#include "gzstream/gzstream.h"

#ifndef DEPLOID_SRC_VCFREADER_HPP_
#define DEPLOID_SRC_VCFREADER_HPP_

// using namespace std;
// using std::endl;

struct InvalidVcf : public InvalidInput{
    explicit InvalidVcf(string str):InvalidInput(str) {
    }
    virtual ~InvalidVcf() throw() {}
    // virtual const char* what () const noexcept {
        // return throwMsg.c_str();
    // }
};


struct VcfInvalidHeaderFieldNames : public InvalidVcf{
    explicit VcfInvalidHeaderFieldNames(string str1, string str2):
                InvalidVcf(str1) {
        this->reason = " VCF field header expects: ";
        throwMsg = this->reason + this->src + ", " + str2 + " was found!";
    }
    ~VcfInvalidHeaderFieldNames() throw() {}
};


struct VcfInvalidVariantEntry : public InvalidVcf{
    explicit VcfInvalidVariantEntry(string str):InvalidVcf(str) {}
    virtual ~VcfInvalidVariantEntry() throw() {}
    // virtual const char* what () const noexcept {
        // return throwMsg.c_str();
    // }
};


struct VcfCoverageFieldNotFound : public VcfInvalidVariantEntry{
    explicit VcfCoverageFieldNotFound(string str):VcfInvalidVariantEntry(str) {
        this->reason = "Coverage field AD was not found in the FORMAT, found: ";
        throwMsg = this->reason + this->src;
    }
    ~VcfCoverageFieldNotFound() throw() {}
};


struct VcfVQSLODNotFound : public VcfInvalidVariantEntry{
    explicit VcfVQSLODNotFound(string str):VcfInvalidVariantEntry(str) {
        this->reason = "VQSLOD was note found, check: ";
        throwMsg = this->reason + this->src;
    }
    ~VcfVQSLODNotFound() throw() {}
};


class VariantLine{
  friend class VcfReader;
  friend class DEploidIO;
 public:
    explicit VariantLine(string tmpLine, size_t sampleColumnIndex, bool extractPlaf = false);
    ~VariantLine() {}

 private:
    string tmpLine_;
    string tmpStr_;

    void init(string tmpLine, size_t sampleColumnIndex, bool extractPlaf);

    void extract_field_CHROM();
    void extract_field_POS();
    void extract_field_ID();
    void extract_field_REF();
    void extract_field_ALT();
    void extract_field_QUAL();
    void extract_field_FILTER();
    void extract_field_INFO();
    void extract_field_FORMAT();
    void extract_field_VARIANT();

    size_t feildStart_;
    size_t fieldEnd_;
    size_t fieldIndex_;

    string chromStr;
    string posStr;
    string idStr;
    string refStr;
    string altStr;
    string qualStr;
    string filterStr;
    string infoStr;
    string formatStr;
    int adFieldIndex_;

    int ref;
    int alt;
    double vqslod;
    double plaf;
    size_t sampleColumnIndex_;
    bool extractPlaf_;
};



/*! \brief VCF file reader @ingroup group_data */
class VcfReader : public VariantIndex {
#ifdef UNITTEST
  friend class TestVCF;
#endif
  friend class DEploidIO;
 public:
    // Constructors and Destructors
    explicit VcfReader(string fileName, string sampleName, bool extractPlaf = false);
    // parse in exclude sites
    ~VcfReader() {}

    // Members and Methods
    vector <string> headerLines;  // calling from python, need to be public
    vector <double> refCount;  // calling from python, need to be public
    vector <double> altCount;  // calling from python, need to be public
    vector <double> vqslod;  // calling from python, need to be public
    vector <double> plaf;
    void finalize();  // calling from python, need to be public

 private:
    vector <VariantLine> variants;
    vector <VariantLine> keptVariants;
    vector <size_t> legitVqslodAt;
    string fileName_;
    ifstream inFile;
    igzstream inFileGz;
    bool isCompressed_;
    bool isCompressed() const { return this->isCompressed_; }
    void setIsCompressed(const bool compressed) {
        this->isCompressed_ = compressed; }
    void checkFileCompressed();
    string sampleName_;
    size_t sampleColumnIndex_;
    string tmpLine_;
    string tmpStr_;
    bool extractPlaf_;

    // Methods
    void init(string fileName);
    void readVariants();
    void readHeader();
    void checkFeilds();
    void findLegitSnpsGivenVQSLOD(double vqslodThreshold);
    void findLegitSnpsGivenVQSLODHalf(double vqslodThreshold);

    void getChromList();
    void removeMarkers();

    // Debug tools
    bool printSampleName();
};

#endif  // DEPLOID_SRC_VCFREADER_HPP_
