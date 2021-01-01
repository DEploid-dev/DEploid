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


#include <iomanip>        // std::setw
#include "dEploidIO.hpp"


void DEploidIO::operation_printVersion(std::ostream& out) {
    out << endl
        << "dEploid " << VERSION
        << endl
        << "Git commit (DEploid): " << dEploidGitVersion_ << endl
        << "Git commit (Lasso): " << lassoGitVersion_ << endl;
}


void DEploidIO::operation_printHelp(std::ostream& out) {
    out << endl
        << "dEploid " << VERSION << endl
        << endl;
    out << "Contact: Joe Zhu <sha.joe.zhu@gmail.com>" << endl
        << endl;
    out << "Usage:"
        << endl;
    out << setw(20) << "-h or -help"         << "  --  "
        << "Help. List the following content." << endl;
    out << setw(20) << "-v or -version"      << "  --  "
        << "DEploid version." << endl;
    out << setw(20) << "-vcf STR"            << "  --  "
        << "VCF file path." << endl;
    out << setw(20) << "-ref STR"            << "  --  "
        << "File path of reference allele count." << endl;
    out << setw(20) << "-alt STR"            << "  --  "
        << "File path of alternative allele count." << endl;
    out << setw(20) << "-plaf STR"           << "  --  "
        << "File path of population level allele frequencies." << endl;
    out << setw(20) << "-panel STR"          << "  --  "
        << "File path of the reference panel." << endl;
    out << setw(20) << "-exclude STR"        << "  --  "
        << "File path of sites to be excluded." << endl;
    out << setw(20) << "-o STR"              << "  --  "
        << "Specify the file name prefix of the output." << endl;
    out << setw(20) << "-p INT"              << "  --  "
        << "Out put precision (default value 8)." << endl;
    out << setw(20) << "-k INT"              << "  --  "
        << "Number of strain (default value 5)." << endl;
    out << setw(20) << "-seed INT"           << "  --  "
        << "Random seed." << endl;
    out << setw(20) << "-nSample INT"        << "  --  "
        << "Number of MCMC samples." << endl;
    out << setw(20) << "-rate INT"           << "  --  "
        << "MCMC sample rate." << endl;
    out << setw(20) << "-noPanel"            << "  --  "
        << "Use population level allele frequency as prior." << endl;
    out << setw(20) << "-forbidUpdateProp"   << "  --  "
        << "Forbid MCMC moves to update proportions." << endl;
    out << setw(20) << "-forbidUpdateSingle" << "  --  "
        << "Forbid MCMC moves to update single haplotype." << endl;
    out << setw(20) << "-forbidUpdatePair"   << "  --  "
        << "Forbid MCMC moves to update pair haplotypes." << endl;
    out << setw(20) << "-initialP FLT ..."   << "  --  "
        << "Initialize proportions." << endl;
    out << setw(20) << "-ibd"                << "  --  "
        << "Use DEploid-IBD" << endl;
    out << setw(20) << "-lasso"              << "  --  "
        << "Use DEploid-LASSO" << endl;
    out << setw(20) << "-best"               << "  --  "
        << "Use DEploid-Best-practice" << endl;
    out << endl;
    out << "Note: Please `man docs/_build/man/dEploid.1' for the manual."
        << endl;
    out << endl;
    out << "Examples:" << endl;
    out << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf "
        << "-plaf data/testData/labStrains.test.PLAF.txt "
        << "-o PG0390-CNopanel -noPanel" << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf "
        << "-exclude data/testData/labStrains.test.exclude.txt "
        << "-plaf data/testData/labStrains.test.PLAF.txt "
        << "-o PG0390-CNopanelExclude -noPanel" << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf "
        << "-exclude data/testData/labStrains.test.exclude.txt "
        << "-plaf data/testData/labStrains.test.PLAF.txt "
        << "-o PG0390-C-ibd -panel data/testData/labStrains.test.panel.txt "
        << "-ibd" << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf "
        << "-exclude data/testData/labStrains.test.exclude.txt "
        << "-plaf data/testData/labStrains.test.PLAF.txt "
        << "-o PG0390-C-best -panel data/testData/labStrains.test.panel.txt "
        << "-best" << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf.gz "
        << "-exclude data/testData/labStrains.test.exclude.txt.gz "
        << "-plaf data/testData/labStrains.test.PLAF.txt.gz -o PG0390-C-best "
        << "-panel data/testData/labStrains.test.panel.txt.gz -best" << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf "
        << "-exclude data/testData/labStrains.test.exclude.txt "
        << "-plaf data/testData/labStrains.test.PLAF.txt "
        << "-o PG0390-CPanelExclude "
        << "-panel data/testData/labStrains.test.panel.txt "
        << "-painting PG0390-CPanelExclude.hap" << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf "
        << "-plaf data/testData/labStrains.test.PLAF.txt "
        << "-o PG0390-CNopanel -noPanel -k 2 -nSample 250 -rate 8 -burn 0.67 "
        << "-ibd "  << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf "
        << "-plaf data/testData/labStrains.test.PLAF.txt "
        << "-o PG0390-CNopanel -initialP 0.2 0.8 -ibdPainting"  << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf "
        << "-exclude data/testData/labStrains.test.exclude.txt "
        << "-plaf data/testData/labStrains.test.PLAF.txt "
        << "-o PG0390-CLassoExclude "
        << "-panel data/testData/labStrains.test.panel.txt "
        << "-initialP 0.2 0.8 -writePanel -lasso" << endl;
    out << "./dEploid -vcf data/testData/PG0390-C.test.vcf.gz "
        << "-sample PG0390-C -plafFromVcf "
        << "-exclude data/testData/labStrains.test.exclude.txt.gz "
        << "-o PG0390-C-best "
        << "-panel data/testData/labStrains.test.panel.txt.gz -best" << endl;
}

