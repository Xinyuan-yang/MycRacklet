#!/bin/bash

svn diff > $1/release.patch
svn --verbose status > $1/full-revision-detail.txt
echo "#include <string>" 
echo "std::string lm_release_info =  \" \\n \\" 
svn info  | sed s/\"/\\\\\"/g | sed s/\\\*/\\\\\*/g | sed 's/$/ \\n\\/'
#svn diff  |  sed "s/\\\\/\\\\\\\\/" | sed s/\"/\\\\\"/g  | sed 's/$/ \\n \\/' | sed "s/\\\\\*/\\\\\\\\*/g" | sed "s/\\\\\//\\\\\\\\\//g" 
echo "\";"