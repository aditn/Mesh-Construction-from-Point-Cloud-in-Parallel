/* needed for file processing */
#include <iostream>
#include <fstream>

/* needed for malloc and printf, scanf, cmd line args and what not */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* needed for openMP! */
#include <omp.h>

/* needed for vectors */
#include <vector>

/* our file dependencies */
#include "structures.h"
#include "approximateMesh.h"
#include "parseOBJ.h"
#include "ourTimer.h"
#include "constants.h"

/* linear algebra library */
#include <Eigen/Dense>

using namespace Eigen;
bool DEBUG = false;
int numThreads=1;
float rho=-1;

int main(int argc, char* argv[]){
  char* in_filename = NULL;
  int opt;
#ifdef USE_OMP
  while((opt = getopt(argc,argv,"f:t:dp:")) != -1){
#else
  while((opt = getopt(argc,argv,"f:dp:")) != -1){
#endif
    switch(opt){
      case 'f': //user must specify an input argument
        in_filename = optarg;
        break;
      case 'p':
        rho = atof(optarg);
        printf("using rho of %.3f\n",rho);
        break;
#ifdef USE_OMP
      case 't':
        numThreads = atoi(optarg);
        printf("using %d threads\n",numThreads);
        if(numThreads<1){
          printf("You want %d threads? No, no you don't\n",numThreads);
          return 1;
        }
        break;
#endif
      case 'd':
        DEBUG = true;
        printf("debug flag set!\n");
        break;
      default:
        fprintf(stderr,"Error. Usage is <exec> -f <filename>\n");
#ifdef USE_OMP
        fprintf(stderr,"You can also use -d to debug and -t <int> to set numthreads\n");
#else
        fprintf(stderr,"You can also use -d to debug\n");
#endif
        return 1;
    }
  }
  if(rho<0){
    fprintf(stderr,"You need to set a positive rho value\n");
    return 1;
  }
  omp_set_num_threads(numThreads);
  startTimer();

  //step 1: parse input file/stream

  std::vector<Vector3f> vertices = parseFilePoints(in_filename);

  int numPoints = vertices.size();
  Vector3f* points = (Vector3f*) malloc(sizeof(Vector3f)*numPoints);
  for(int i=0;i<numPoints;i++) points[i] = vertices[i];
  printf("time %.4fs\n",timeSince());
  
  //step 2: create mesh in parallel
  std::vector<Vector3f> finalVertices;
  std::vector<Edge> edges;
  //will pass outputs into function, which will write to them
  approximateMesh(points,numPoints,finalVertices,edges); 
  printf("time %.4fs\n",timeSince());

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
  printf("time %.4fs\n",timeSince());
  return 0;
}
