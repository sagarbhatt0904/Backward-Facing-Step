#!/bin/sh

g++ -O2 -Wall  ./src/main.cpp -lblas -fopenmp -o NS

./NS

rm NS