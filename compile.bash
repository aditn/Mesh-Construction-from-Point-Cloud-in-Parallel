#!/bin/bash
g++ -Wall -I util/ main.cpp tangentPlane.cpp parseOBJ.cpp approximateMesh.cpp -o runnable -fopenmp
