/*I mean this could be c or c++ but I feel like c++ is nicer and lets us use cooler data structures*/

/* needed for file processing */
#include <iostream>
#include <fstream>

/* needed for vectors */
#include <vector>

/* needed for malloc and printf and scanf and what not */
#include <stdio.h>
#include <stdlib.h>

/* sqrt and square */
#include <cmath>

/* data structures */
#include "structures.h"
#include "tangentPlane.h"

/* parseFile takes a filename as input, tests it as an obj, and grabs all of its vertices
 * It returns of Point vector of these vertices
 * An error will be thrown if file cannot be read */
std::vector<V3> parseFile(const char* filename){
  std::vector<V3> points;
  std::ifstream fin;
  fin.open(filename);
  if(!fin.good()){
    printf("Couldn't open input file!\n");
    exit(0);
  }else{
    char line[MAX_LINE_SIZE],line_type[MAX_LINE_SIZE];
    while(fin.getline(line,MAX_LINE_SIZE)){
      float x,y,z;
      sscanf(line,"%s %f %f %f",line_type,&x,&y,&z);
      printf("line of type %s is (%.3f,%.3f,%.3f)\n",line_type,x,y,z);
      if(line_type[0]=='v'){ //vertex
        V3 p;
        p.x=x;p.y=y;p.z=z;
        points.push_back(p);
      }else printf("couldn't recognize type %s\n",line_type);
    }
  }
  fin.close();
  return points;
}

void saveMesh(V3* points,int numPoints, char* out_filename){
  std::ofstream fout;
  fout.open(out_filename);
  if(!fout.good()){
    printf("couldn't write to %s\n",out_filename);
    exit(0);
  }else{
    for(int i=0;i<numPoints;i++){
      float x,y,z;
      x=points[i].x;
      y=points[i].y;
      z=points[i].z;
      fout << "v " << x << " " << y << " " << z << "\n";
    }
  }
  fout.close();
}

inline void printPoint(V3 p){
  printf("(%.3f,%.3f,%.3f)\n",p.x,p.y,p.z);
}

int main(){
  //for now we can read/write from obj file to keep things simple
  //until we decide how we're actually gonna grab kinect stuff

  //step 1: parse input file/stream
  const char* in_filename = "inputs/test.obj";
  std::vector<V3> vertices = parseFile(in_filename);
  int numPoints = vertices.size();
  V3* points = (V3*) malloc(sizeof(V3)*numPoints);
  for(int i=0;i<numPoints;i++) points[i] = vertices[i];
  for(int i=0;i<numPoints;i++) printPoint(points[i]);
  //step 2: create mesh in parallel
  //step 2a: get neighborhood for each point
  //step 2b: get centroid & PCA normals based off of neighborhoods
  Plane* planes = computeTangentPlanes(points,numPoints,0.6,0.5);
  for(int i=0;i<numPoints;i++){
    printf("Plane %d: \n\t",i);
    printPoint(planes[i].center);
    printf("\t");printPoint(planes[i].normal);
  }
  //step 2c: Propogate normal directions for surface consistency
  //step 2d: Create cubes around these centers
  //step 2e: Approximate mesh based on differences between cubes
  //step 2f: Refine mesh

  //step 3: write mesh to output file
  int slen = 0;
  while(in_filename[slen++]!='\0');
  char out_filename[slen+5];
  for(int i=0;i<3;i++) out_filename[i] = "out"[i];
  for(int i=2;i<slen-5;i++) out_filename[i+1]=in_filename[i];
  for(int i=0;i<9;i++) out_filename[slen-4+i]="_out.obj"[i]; //9 includes null terminator
  printf("saving to %s\n",out_filename);
  saveMesh(points,numPoints,out_filename);
  return 0;
}
