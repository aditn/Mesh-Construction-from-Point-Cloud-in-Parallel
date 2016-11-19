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
  float dot(V3 v){
    return this->x*v.x+this->y*v.y+this->z*v.z;
  }
  void scale(float val){
    this->x*=val;
    this->y*=val;
    this->z*=val;
  }
  float dist(V3 v){
    return std::pow(std::pow(this->x-v.x,2)+std::pow(this->y-v.y,2)+std::pow(this->z-v.z,2),0.5);
  } 
  V3 cross(V3 v){
    V3 out;
    out.x = this->y*v.z-this->z*v.y;
    out.y = this->z*v.x-this->x*v.z;
    out.z = this->x*v.y-this->y*v.x;
    return out;
  }
};
inline V3 operator-(const V3& a, const V3& b){
    V3 out;
    out.x = a.x-b.x;
    out.y = a.y-b.y;
    out.z = a.z-b.z;
    return out;
}

typedef struct{
  V3 center, normal;
} Plane;
