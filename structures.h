#include <cmath>
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

typedef struct{
  V3 center, normal;
} Plane;
