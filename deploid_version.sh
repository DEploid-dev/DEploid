#!/bin/sh
git show HEAD | head -1 | sed -e "s/commit //g" | cat
