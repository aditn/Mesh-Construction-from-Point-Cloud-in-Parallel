#include <cmath>
#include <stdio.h>
#define MAX_LINE_SIZE 1024 //our max line size for obj file

struct V3{
  float x,y,z;
  V3(){
    this->x=0;
    this->y=0;
    this->z=0;
  }
  V3(float x,float y,float z){
    this->x = x;
    this->y = y;
    this->z = z;
  }
  void add(V3 v){
    this->x+=v.x;
    this->y+=v.y;
    this->z+=v.z;
  }
  V3 dot(V3 v){
    return V3(this->x*v.x,this->y*v.y,this->z*v.z);
  }
  void scale(float val){
    this->x*=val;
    this->y*=val;
    this->z*=val;
  }
  float dist(V3 v){
    return std::pow(std::pow(this->x-v.x,2)+std::pow(this->y-v.y,2)+std::pow(this->z-v.z,2),0.5);
  }
};

struct Plane{
  V3 center, normal;
};

struct mat{
  short row, col;
  double** matrix;
  mat(){
    this->row=0;
    this->col=0;
  }
  mat(int row, int col, double** pointerToMat){
    this->row = row;
    this->col = col;
    this->matrix = pointerToMat;
  }
};
