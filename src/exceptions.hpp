/*
 * pfDeconv is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016, Sha (Joe) Zhu, Jacob Almagro and Prof. Gil McVean
 *
 * This file is part of pfDeconv.
 *
 * scrm is free software: you can redistribute it and/or modify
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

*/


#include <exception>

#ifndef EXCEPTION
#define EXCEPTION

using namespace std;

struct ooops : std::exception {
  const char* what() const noexcept {return "Ooops!\n";}
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
    this->src      = str;
    this->reason   = "";
  }
  virtual ~InvalidInput() {}
  virtual const char* what () const noexcept {
    return throwMsg.c_str();
  }
};


struct NotEnoughArg : public InvalidInput{
  NotEnoughArg( string str ):InvalidInput( str ){
    this->reason = "Not enough parameters when parsing option: ";
    throwMsg = this->reason + this->src;
  }
  ~NotEnoughArg(){}
};


struct WrongType : public InvalidInput{
  WrongType( string str ):InvalidInput( str ){
    this->reason = "Wrong type for parsing: ";
    throwMsg = this->reason + this->src;
  }
  ~WrongType(){}
};


struct InvalidInputFile : public InvalidInput{
  InvalidInputFile( string str ):InvalidInput( str ){
    this->reason = "Invalid input file: ";
    throwMsg = this->reason + this->src;
  }
  ~InvalidInputFile(){}
};


struct FileNameMissing : public InvalidInput{
  FileNameMissing( string str ):InvalidInput( str ){
    this->reason = " file path missing!";
    throwMsg = this->src + this->reason ;
  }
  ~FileNameMissing(){}
};


struct UnknowArg : public InvalidInput{
  UnknowArg( string str ):InvalidInput( str ){
    this->reason = "Unknow option: ";
    throwMsg = this->reason + this->src;
  }
  ~UnknowArg(){}
};


struct FlagsConflict : public InvalidInput{
  FlagsConflict( string str1, string str2 ):InvalidInput( str1 ){
    this->reason = "Flag: ";
    throwMsg = this->reason + this->src + string(" conflict with flag ") + str2;
  }
  ~FlagsConflict(){}
};


struct LociNumberUnequal : public InvalidInput{
  LociNumberUnequal( string str ):InvalidInput( str ){
    this->reason = "Number of sites was wrong (compared to ref count) in: ";
    throwMsg = this->reason + this->src ;
  }
  ~LociNumberUnequal(){}
};


#endif
