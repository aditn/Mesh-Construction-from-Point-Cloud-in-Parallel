/*I mean this could be c or c++ but I feel like c++ is nicer and lets us use cooler data structures*/

#include <iostream>
#include <fstream>

#include <vector>
#include <stdio.h>
#include <stdlib.h>

#define MAX_LINE_SIZE 1024

typedef struct{
  float x,y,z;
} Point;

std::vector<Point> parseFile(const char* filename){
  std::vector<Point> points;
  std::ifstream fin;
  fin.open(filename);
  if(!fin.good()){
    printf("Couldn't open input file!\n");
    exit(0);
  }else{
    while(!fin.eof()){
      char line[MAX_LINE_SIZE],line_type[MAX_LINE_SIZE];
      float x,y,z;
      fin.getline(line,MAX_LINE_SIZE);
      sscanf(line,"%s %f %f %f",line_type,&x,&y,&z);
      printf("line of type %s is (%.3f,%.3f,%.3f)\n",line_type,x,y,z);
      if(line_type[0]=='v'){ //vertex
        Point p;
        p.x=x;p.y=y;p.z=z;
        points.push_back(p);
      }else printf("couldn't recognize type %s\n",line_type);
    }
  }
  fin.close();
  return points;
}

void printPoint(Point p){
  printf("(%.3f,%.3f,%.3f)\n",p.x,p.y,p.z);
}

int main(){
  //for now we can read/write from obj file to keep things simple
  //until we decide how we're actually gonna grab kinect stuff

  //step 1: parse input file/stream
  const char* filename = "test.obj";
  std::vector<Point> vertices = parseFile(filename);
  for(int i=0;i<(int)vertices.size();i++) printPoint(vertices[i]);
  //step 2: create mesh in parallel
  //step 3: write mesh to output file

  return 0;
}
