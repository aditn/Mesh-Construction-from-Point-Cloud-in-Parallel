#!/bin/bash
g++ -Wall -std=c++11 -I util/ main.cpp tangentPlane.cpp parseOBJ.cpp approximateMesh.cpp -o runnable -fopenmp
