#!/bin/sh
git submodule status | grep DEploid-Lasso-lib | sed -e "s/ //g" -e "s/DEploid.*//"
