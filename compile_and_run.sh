#!/bin/bash
cd bin
touch NULL
rm -rf main *.h5 
cd ../build
rm -rf CMakeCache.txt  CMakeFiles  cmake_install.cmake 
cmake ..
make
cd ..
cd bin
echo "-----------------"
./main
cd ..