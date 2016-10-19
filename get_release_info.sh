#!/bin/bash

echo "#include <string>" 
echo "std::string cR_release_info =  \" \\n \\" 
echo "Remote origin info: \\n \\"
git config --get remote.origin.url | sed 's/\\/bed/g' | sed 's/$/ \\n\\/'
git log origin/master -1 --pretty=fuller  | sed 's/\\/bed/g' | sed 's/$/ \\n\\/'
echo "\\n \\"
echo "Local branch info: \\n \\"
git log -1 --pretty=fuller  | sed 's/\\/bed/g' | sed 's/$/ \\n\\/'    
echo "\\n \\"
echo "Local branch status: \\n \\"
git status -uno  | sed 's/\\/bed/g' |  sed "s/\\\\/\\\\\\\\/" | sed s/\"/\\\\\"/g  | sed 's/$/ \\n \\/' | sed "s/\\\\\*/\\\\\\\\*/g" | sed "s/\\\\\//\\\\\\\\\//g" 
echo "\\n \\"
echo "\\n \\"
echo "-> More details about the version of cRacklet sources can be found in file Sources_info.cra. \\n \\"
echo "\";"
echo "#include <string>" 
echo "std::string cR_sources_info =  \" \\n \\" 
git diff | sed 's/\\/bed/g' |  sed "s/\\\\/\\\\\\\\/" | sed s/\"/\\\\\"/g  | sed 's/$/ \\n \\/' | sed "s/\\\\\*/\\\\\\\\*/g" | sed "s/\\\\\//\\\\\\\\\//g" 
echo "\";"
