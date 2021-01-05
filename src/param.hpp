/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2020 - 2021  Sha (Joe) Zhu
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


#ifndef DEPLOID_SRC_PARAM_HPP_
#define DEPLOID_SRC_PARAM_HPP_

template <class T> class Parameter {
friend class DEploidIO;

 private:
    bool useDefault_;
    void setUseDefault(const bool setTo) { this->useDefault_ = setTo;}
    bool useDefault() const { return this->useDefault_; }

    bool useBest_;
    void setUseBest(const bool setTo) { this->useBest_ = setTo;}
    bool useBest() const { return this->useBest_; }

    bool useUserDefined_;
    void setUseUserDefined(const bool setTo) {this->useUserDefined_ = setTo;}
    bool useUserDefined() const { return this->useUserDefined_; }

    T default_;
    void setDefault(const T setTo) {this->default_ = setTo;}
    T best_;
    void setBest(const T setTo) {
        this->setUseBest(true);
        this->best_ = setTo;
    }

    T userDefined_;
    void setUserDefined(const T setTo) {
        this->setUseUserDefined(true);
        this->userDefined_ = setTo;
    }

    void init(T value) {
        this->setDefault(value);
        this->setBest(value);
        this->setUserDefined(value);
        this->setUseDefault(true);
        this->setUseBest(false);
        this->setUseUserDefined(false);
    }

    void makeCopy(const Parameter <T> &currentParam) {
        this->setDefault(currentParam.default_);
        this->setBest(currentParam.best_);
        this->setUserDefined(currentParam.userDefined_);
        // Make sure to set the value states after set values
        this->setUseDefault(currentParam.useDefault());
        this->setUseBest(currentParam.useBest());
        this->setUseUserDefined(currentParam.useUserDefined());
    }

 public:
    Parameter <T> () { }
    Parameter <T> (T value) {
        this->init(value);
    }
    T getValue();
    Parameter <T> (const Parameter <T> &currentParam) {
        this->makeCopy(currentParam);
    }
    ~Parameter() { }
};

// Parameter <double> a_Parameter_class;

template <class T> T Parameter<T>::getValue() {
    if (useUserDefined_) {
        return userDefined_;
    } else {
        if (useBest_) {
            return best_;
        } else {
            // assert(useDefault_);
            return default_;
        }
    }
}

#endif  // DEPLOID_SRC_PARAM_HPP_

