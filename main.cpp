/* needed for file processing */
#include <iostream>
#include <fstream>

/* needed for malloc and printf, scanf, cmd line args and what not */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <chrono> //for timing

/* needed for openMP! */
#include <omp.h>

/* needed for vectors */
#include <vector>

/* our file dependencies */
#include "structures.h"
#include "approximateMesh.h"
#include "parseOBJ.h"

#define rho 0.4 //given rho radius sphere, at least one point?
#define delta 0.6 //noise
bool DEBUG = false;
int numThreads;


inline float timeSince(std::chrono::time_point<std::chrono::steady_clock> start){
  float ms = std::chrono::duration<float,std::milli>(std::chrono::steady_clock::now()-start).count();
  return ms*0.001;
}

int main(int argc, char* argv[]){
  char* in_filename = NULL;
  int opt;
  numThreads=1;
  while((opt = getopt(argc,argv,"f:t:d")) != -1){
    switch(opt){
      case 'f': //user must specify an input argument
        in_filename = optarg;
        break;
      case 't':
        numThreads = atoi(optarg);
        printf("using %d threads\n",numThreads);
        if(numThreads<1){
          printf("You want %d threads? No, no you don't\n",numThreads);
          return 1;
        }
        break;
      case 'd':
        DEBUG = true;
        printf("debug flag set!\n");
        break;
      default:
        fprintf(stderr,"Error. Usage is <exec> -f <filename>\n");
        fprintf(stderr,"You can also use -d to debug and -t <int> to set numthreads\n");
        return 1;
    }
  }
  omp_set_num_threads(numThreads);
  auto start = std::chrono::steady_clock::now();

  //step 1: parse input file/stream
  std::vector<V3> vertices = parseFilePoints(in_filename);
  int numPoints = vertices.size();
  V3* points = (V3*) malloc(sizeof(V3)*numPoints);
  for(int i=0;i<numPoints;i++) points[i] = vertices[i];
  printf("time %.4fs\n",timeSince(start));
  
  //step 2: create mesh in parallel
  std::vector<V3> finalVertices;
  std::vector<Edge> edges;
  //will pass outputs into function, which will write to them
  approximateMesh(points,numPoints,rho,delta,finalVertices,edges); 
  printf("time %.4fs\n",timeSince(start));

  //step 2f: Refine mesh

  //step 3: write mesh to output file
  int slen = 0;
  while(in_filename[slen++]!='\0');//get length of string
  char out_filename[slen+5];
  for(int i=0;i<3;i++) out_filename[i] = "out"[i];//change inputs/X to outputs/X
  for(int i=2;i<slen-5;i++) out_filename[i+1]=in_filename[i];
  for(int i=0;i<9;i++) out_filename[slen-4+i]="_out.obj"[i]; //9 includes null terminator
  printf("saving to %s\n",out_filename);
  saveMesh(finalVertices,edges,out_filename);
  printf("time %.4fs\n",timeSince(start));
  return 0;
}
