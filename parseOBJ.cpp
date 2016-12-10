/* file processing */
#include <iostream>
#include <fstream>

/* malloc,printf,scanf,... */
#include <stdio.h>
#include <stdlib.h>

/* vectors */
#include <vector>

/* our datatypes */
#include "structures.h"
#include "parseOBJ.h"

/* Eigen Library */
#include <Eigen/Dense>

using namespace Eigen;

#define MAX_LINE_SIZE 1024 //our max line size for obj file

std::vector<Vector3f> parseFilePoints(const char* filename){
  std::vector<Vector3f> points;
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
      if(line_type[0]=='v') points.push_back(Vector3f(x,y,z)); //vertices are points in point cloud
      else printf("couldn't recognize type %s\n",line_type);
    }
  }
  fin.close();
  return points;
}

void saveMesh(std::vector<Vector3f> V, std::vector<Edge> edges, const char* out_filename){
  std::ofstream fout;
  fout.open(out_filename);
  if(!fout.good()){
    printf("couldn't write to %s\n",out_filename);
    exit(0);
  }else{
    //print all of the vertices
    for(unsigned int i=0;i<V.size();i++){
      fout << "v " << V[i](0) << " " << V[i](1) << " " << V[i](2) << "\n";
    }

    //obj files can have "l" datatype as a line
    //this is how we will represent edges
    for(unsigned int i=0;i<edges.size();i++){
      fout << "l " << (edges[i].v1+1) << " " << (edges[i].v2+1) << "\n";
    }
  }
  fout.close();
}
