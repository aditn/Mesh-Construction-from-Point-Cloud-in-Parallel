/* needed for file processing */
#include <iostream>
#include <fstream>

/* needed for malloc and printf, scanf, cmd line args and what not */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* needed for vectors */
#include <vector>

/* our file dependencies */
#include "structures.h"
#include "approximateMesh.h"

#define rho 0.4
#define delta 0.6

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
      if(line_type[0]=='v') points.push_back(V3(x,y,z)); //vertices are points in point cloud
      else printf("couldn't recognize type %s\n",line_type);
    }
  }
  fin.close();
  return points;
}

void saveMesh(std::vector<V3> V, std::vector<Edge> edges, char* out_filename){
  std::ofstream fout;
  fout.open(out_filename);
  if(!fout.good()){
    printf("couldn't write to %s\n",out_filename);
    exit(0);
  }else{
    //print all of the vertices
    for(unsigned int i=0;i<V.size();i++){
      fout << "v " << V[i].x << " " << V[i].y << " " << V[i].z << "\n";
    }

    //obj files can have "l" datatype as a line
    //this is how we will represent edges
    for(unsigned int i=0;i<edges.size();i++){
      fout << "l " << (edges[i].v1+1) << " " << (edges[i].v2+1) << "\n";
    }
  }
  fout.close();
}

int main(int argc, char* argv[]){
  char* in_filename = NULL;
  int opt;
  while((opt = getopt(argc,argv,"f:")) != -1){
    switch(opt){
      case 'f': //user must specify an input argument
        in_filename = optarg;
        break;
      default:
        fprintf(stderr,"Error. Usage is <exec> -f <filename>\n");
        return 1;
    }
  }

  //step 1: parse input file/stream
  std::vector<V3> vertices = parseFile(in_filename);
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
