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

#define rho 0.4 //given rho radius sphere, at least one point?
#define delta 0.6 //noise
bool DEBUG = false;

int main(int argc, char* argv[]){
  char* in_filename = NULL;
  int opt;
  while((opt = getopt(argc,argv,"f:d")) != -1){
    switch(opt){
      case 'f': //user must specify an input argument
        in_filename = optarg;
        break;
      case 'd':
        DEBUG = true;
        printf("debug flag set!\n");
        break;
      default:
        fprintf(stderr,"Error. Usage is <exec> -f <filename>\n");
        return 1;
    }
  }
  omp_set_num_threads(10);
  #pragma omp parallel for
  for(int i=0;i<10;i++) printf("hello from %d\n",i);

  //step 1: parse input file/stream
  std::vector<V3> vertices = parseFilePoints(in_filename);
  int numPoints = vertices.size();
  V3* points = (V3*) malloc(sizeof(V3)*numPoints);
  for(int i=0;i<numPoints;i++) points[i] = vertices[i];
  
  //step 2: create mesh in parallel
  std::vector<V3> finalVertices;
  std::vector<Edge> edges;
  //will pass outputs into function, which will write to them
  approximateMesh(points,numPoints,rho,delta,finalVertices,edges); 

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
  return 0;
}
