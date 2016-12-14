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
DebugInfo debug_info;

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
        break;
#ifdef USE_OMP
      case 't':
        numThreads = atoi(optarg);
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
  if(rho<=0){
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
  debug_info.fileParse = timeSince();  

  //step 2: create mesh in parallel
  std::vector<Vector3f> finalVertices;
  std::vector<Edge> edges;
  //will pass outputs into function, which will write to them
  approximateMesh(points,numPoints,finalVertices,edges); 

  //step 2f: Refine mesh

  //step 3: write mesh to output file
  int slen = 0;
  while(in_filename[slen++]!='\0');//get length of string
  char out_filename[slen+5];
  for(int i=0;i<3;i++) out_filename[i] = "out"[i];//change inputs/X to outputs/X
  for(int i=2;i<slen-5;i++) out_filename[i+1]=in_filename[i];
  for(int i=0;i<9;i++) out_filename[slen-4+i]="_out.obj"[i]; //9 includes null terminator
  saveMesh(finalVertices,edges,out_filename);
  debug_info.complete = timeSince();

  printf("-------File Stats-----------\n");
  printf("rho is %.3f\n",rho);
  printf("running on %d processors\n",numThreads);
  printf("point cloud bbox is (%.3f,%.3f,%.3f) to (%.3f,%.3f,%.3f)\n",
    debug_info.cloud.min(0),debug_info.cloud.min(1),debug_info.cloud.min(2),
    debug_info.cloud.max(0),debug_info.cloud.max(1),debug_info.cloud.max(2));
  printf("space quantized to %d cubes (%dx%dx%d)\n",
    debug_info.cubex*debug_info.cubey*debug_info.cubez,
    debug_info.cubex,debug_info.cubey,debug_info.cubez);
  printf("%.3f%% redundancy eliminated\n",debug_info.redundancy);
  printf("--------timing info---------\n");
  printf("%.4fs: file parsed\n",debug_info.fileParse);
  printf("%.4fs: planes created\n",debug_info.planeCreation);
  printf("%.4fs: %d-Nearest Neighbor graph formed\n",debug_info.kNN,NUM_NEIGHBORS);
  printf("%.4fs: MST Generated and normals propagated\n",debug_info.MST);
  printf("%.4fs: cube corner distances computed\n",debug_info.distsFound);
  printf("%.4fs: mesh points and edges computed\n",debug_info.mesh);
  printf("%.4fs: mesh points and edges deduped\n",debug_info.dedup);
  printf("%.4fs: mesh written to file.\n",debug_info.complete);
  printf("----------------------------\n");
  printf("Result saved to %s\n",out_filename);
  return 0;
}
