#!/bin/bash

git diff > $1/release.patch
git show > $1/full-revision-detail.txt
echo "#include <string>" 
echo "std::string lm_release_info =  \" \\n \\" 
echo "Local branch status: \\n \\"
git log -1 --pretty=fuller  | sed s/\"/\\\\\"/g | sed s/\\\*/\\\\\*/g | sed 's/$/ \\n\\/'
echo "\\n \\"    
echo "Remote origin status: \\n \\"
git config --get remote.origin.url | sed s/\"/\\\\\"/g | sed s/\\\*/\\\\\*/g | sed 's/$/ \\n\\/'
git log origin/master -1 --pretty=fuller  | sed s/\"/\\\\\"/g | sed s/\\\*/\\\\\*/g | sed 's/$/ \\n\\/'
echo "\\n \\"
git status -uno  |  sed "s/\\\\/\\\\\\\\/" | sed s/\"/\\\\\"/g  | sed 's/$/ \\n \\/' | sed "s/\\\\\*/\\\\\\\\*/g" | sed "s/\\\\\//\\\\\\\\\//g" 
echo "\";"
