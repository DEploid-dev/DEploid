#!/bin/bash

#if [ ! -z "$src_dir+x" ]; then
src_dir=$PWD/src
#fi

if [ ! -f "${src_dir}/dEploid.cpp" ]; then
 "Error: please execute \".checkfiles.sh\" from DEploid root directory"
 exit 1
fi


OS=`uname -s`
#echo "operating system is ${OS}"
# VERSION=$(cat VERSION | tr -d '\n')
VERSION=v0.1-release
echo -n "." # echo "ok not git clone, now check plot and figure"

#rm master.tar.gz  # just making sure download the one we actually want to

extract_submodule_from(){
 account=$1
 submodule=$2
 plactto=$3
 echo "Extract submodule $submodule from account $account, place at $3"
 if [ ! -d ${submodule}_dir ]; then mkdir ${submodule}_dir; fi
 echo -n "." # echo "ok, do something" #

 if [[ "${OS}" == "Linux" ]]; then
   wget --no-check-certificate https://github.com/${account}/${submodule}/archive/${VERSION}.tar.gz -o /dev/null
 elif [[ "${OS}" == "Darwin"* ]]; then
   curl -LOk https://github.com/${account}/${submodule}/archive/${VERSION}.tar.gz
 else
   echo "Unknown OS, fail to download package from https://github.com/${account}/${submodule}/archive/${VERSION}.tar.gz"
   echo "Please contact Joe at sha.joe.zhu@gmail.com if assistance is needed"
   exit 1
 fi
ls *gz

 if [ ! -f ${VERSION}.tar.gz ]; then
   echo "Error: Download package from https://github.com/${account}/${submodule}/archive/${VERSION}.tar.gz failed"
   echo "Please contact Joe at sha.joe.zhu@gmail.com if assistance is needed"
   exit 1
 fi
 tar -xf ${VERSION}.tar.gz -C ${submodule}_dir
 if [ ! -d ${placeto} ]; then mkdir ${placeto}; fi
 # ls  ${submodule}_dir/*
 cp -r ${submodule}_dir/* ${placeto}/
 rm -r ${submodule}_dir ${VERSION}.tar.gz
}


ACCOUNT=DEploid-dev
SUBMODLUE=DEploid-Utilities
PLACETO=utilities
tmp_dir=$PWD/utilities
if [ -f ${tmp_dir}/DEploidR.R  -a -f ${tmp_dir}/dEploidTools.r ];then
 echo -n "." # echo "ok, do nothing"
else
 extract_submodule_from $ACCOUNT $SUBMODLUE $PLACETO
fi

ACCOUNT=shajoezhu
SUBMODLUE=DEploid-Lasso-lib
PLACETO=$PWD/src/lasso
tmp_dir=$PWD/src/lasso
if [ -f ${tmp_dir}/main.cpp  -a -f ${tmp_dir}/README.md ];then
 echo -n "." # echo "ok, do nothing"
else
 extract_submodule_from $ACCOUNT $SUBMODLUE $PLACETO
fi
