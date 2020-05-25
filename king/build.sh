#!/bin/bash

export CXX=${PREFIX}/bin/g++

mkdir -p $PREFIX/bin
cp ./KINGcode.tar.gz $PREFIX/bin/.
wget http://people.virginia.edu/~wc9c/KING/KINGcode.tar.gz
tar -xzvf KINGcode.tar.gz

g++ -lm -lz -O2 -fopenmp -o king *.cpp
