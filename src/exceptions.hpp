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


#include <string>  /* string */
#include <exception>

#ifndef EXCEPTION
#define EXCEPTION

using namespace std;


struct ShouldNotBeCalled : std::exception{
  explicit ShouldNotBeCalled(){ }
  virtual ~ShouldNotBeCalled() throw() {}
  virtual const char* what () const noexcept {
      return string("Should not reach here").c_str();
  }
};


struct VirtualFunctionShouldNotBeCalled : public ShouldNotBeCalled{
  VirtualFunctionShouldNotBeCalled():ShouldNotBeCalled(){}
  ~VirtualFunctionShouldNotBeCalled() throw() {}
};


struct InvalidInput : std::exception {
  string src;
  string reason;
  string throwMsg;

  InvalidInput( ){
    this->src      = "";
    this->reason   = "";
  }

  explicit InvalidInput( string str ){
    this->src      = "\033[1;31m" + str + "\033[0m";
    this->reason   = "";
  }
  virtual ~InvalidInput() throw() {}
  virtual const char* what () const noexcept {
      return throwMsg.c_str();
  }
};


struct OutOfVectorSize : std::exception{
  explicit OutOfVectorSize(){ }
  virtual ~OutOfVectorSize() throw() {}
  virtual const char* what () const noexcept {
      return string("Out of vector size!").c_str();
  }
};


struct InvalidK : public InvalidInput{
  InvalidK( ):InvalidInput( ){
    this->reason = "k must be at least 2, when using the flag -ibd.";
    throwMsg = this->reason + this->src;
  }
  ~InvalidK() throw() {}
};


struct NotEnoughArg : public InvalidInput{
  NotEnoughArg( string str ):InvalidInput( str ){
    this->reason = "Not enough parameters when parsing option: ";
    throwMsg = this->reason + this->src;
  }
  ~NotEnoughArg() throw() {}
};


struct VcfOutUnSpecified : public InvalidInput{
  VcfOutUnSpecified( string str ):InvalidInput( str ){
    this->reason = "Missing flag \"-vcfOut\".";
    throwMsg = this->reason + this->src;
  }
  ~VcfOutUnSpecified() throw() {}
};


struct WrongType : public InvalidInput{
  WrongType( string str ):InvalidInput( str ){
    this->reason = "Wrong type for parsing: ";
    throwMsg = this->reason + this->src;
  }
  ~WrongType() throw() {}
};


struct BadConversion : public InvalidInput{
  BadConversion( string str1, string str2 ):InvalidInput( str1 ){
    this->reason = "Bad conversion: ";
    throwMsg = this->reason + this->src + ", int expected. Check input file" + str2;
  }
  ~BadConversion() throw() {}
};


struct InvalidInputFile : public InvalidInput{
  InvalidInputFile( string str ):InvalidInput( str ){
    this->reason = "Invalid input file: ";
    throwMsg = this->reason + this->src;
  }
  ~InvalidInputFile() throw() {}
};


struct FileNameMissing : public InvalidInput{
  FileNameMissing( string str ):InvalidInput( str ){
    this->reason = " file path missing!";
    throwMsg = this->src + this->reason ;
  }
  ~FileNameMissing() throw() {}
};


struct UnknowArg : public InvalidInput{
  UnknowArg( string str ):InvalidInput( str ){
    this->reason = "Unknow option: ";
    throwMsg = this->reason + this->src;
  }
  ~UnknowArg() throw() {}
};


struct FlagsConflict : public InvalidInput{
  FlagsConflict( string str1, string str2 ):InvalidInput( str1 ){
    this->reason = "Flag: ";
    throwMsg = this->reason + this->src + string(" conflict with flag ") + str2;
  }
  ~FlagsConflict() throw() {}
};


struct OutOfRange : public InvalidInput{
  OutOfRange( string str1, string str2 ):InvalidInput( str1 ){
    this->reason = "Flag \"";
    throwMsg = this->reason + this->src + string(" ") + str2 + string("\" out of range [0, 1].");
  }
  ~OutOfRange() throw() {}
};


struct LociNumberUnequal : public InvalidInput{
  LociNumberUnequal( string str ):InvalidInput( str ){
    this->reason = "Number of sites was wrong (compared to ref count) in: ";
    throwMsg = this->reason + this->src ;
  }
  ~LociNumberUnequal() throw() {}
};


struct SumOfPropNotOne : public InvalidInput{
  SumOfPropNotOne( string str ):InvalidInput( str ){
    this->reason = "Sum of initial proportion is not equal to 1, but equals ";
    throwMsg = this->reason + this->src ;
  }
  ~SumOfPropNotOne() throw() {}
};


struct NumOfPropNotMatchNumStrain : public InvalidInput{
  NumOfPropNotMatchNumStrain( string str ):InvalidInput( str ){
    this->reason = "Number of initial proportion do not match number of strains!";
    throwMsg = this->reason + this->src ;
  }
  ~NumOfPropNotMatchNumStrain() throw() {}
};


struct InitialPropUngiven : public InvalidInput{
  InitialPropUngiven( string str ):InvalidInput( str ){
    this->reason = "Initial proportion was not specified.";
    throwMsg = this->reason + this->src ;
  }
  ~InitialPropUngiven() throw() {}
};

#endif
