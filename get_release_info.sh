#!/bin/bash

git diff > $1/release.patch
git show > $1/full-revision-detail.txt
echo "#include <string>" 
echo "std::string lm_release_info =  \" \\n \\" 
git log -1 --pretty=fuller  | sed s/\"/\\\\\"/g | sed s/\\\*/\\\\\*/g | sed 's/$/ \\n\\/'
git status -uno  |  sed "s/\\\\/\\\\\\\\/" | sed s/\"/\\\\\"/g  | sed 's/$/ \\n \\/' | sed "s/\\\\\*/\\\\\\\\*/g" | sed "s/\\\\\//\\\\\\\\\//g" 
echo "\";"
