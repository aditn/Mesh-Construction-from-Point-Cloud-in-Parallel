#include <cmath>

/* linear algebra library */
#include <Eigen/Dense>

#include <stdio.h>
#define TINYNUM 0.000001


/*struct V3{
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
  void sub(V3 v){
    this->x-=v.x;
    this->y-=v.y;
    this->z-=v.z;
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
  V3 unit(){
    float mag = std::pow(this->dot(*this),0.5);
    V3 retval = V3(*this);
    retval.scale(1.0f/mag);
    return retval;
  }
  void print(){
    printf("(%.3f,%.3f,%.3f)\n",this->x,this->y,this->z);
  }
};
inline V3 operator-(const V3 a, const V3 b){
  V3 out;
  out.x = a.x-b.x;
  out.y = a.y-b.y;
  out.z = a.z-b.z;
  return out;
}
inline V3 operator+(const V3 a, const V3 b){
    V3 out;
    out.x = a.x+b.x;
    out.y = a.y+b.y;
    out.z = a.z+b.z;
    return out;
}
inline void operator-=(V3& a, const V3 b){
    a.sub(b);
}
inline void operator+=(V3& a, const V3 b){
    a.add(b);
}
inline V3 operator*(const V3 a,const float b){
  V3 out = V3(a);
  out.scale(b);
  return out;
}
inline V3 operator*(const float b,const V3 a){
  V3 out = V3(a);
  out.scale(b);
  return out;
}
inline void operator*=(V3& a,const float b){
  a.scale(b);
}

inline bool operator==(const V3 a, const V3 b){
  return ((fabs(a.x-b.x)<TINYNUM) && (fabs(a.y-b.y)<TINYNUM)) && (fabs(a.z-b.z)<TINYNUM);
}*/

/* Bounding Box that represents cube */
struct bbox{
  Eigen::Vector3f min,max;
  bbox(){
    this->min = Eigen::Vector3f(0.f,0.f,0.f);
    this->max = Eigen::Vector3f(0.f,0.f,0.f);
  }
  bbox(Eigen::Vector3f p1, Eigen::Vector3f p2){
    this->min = Eigen::Vector3f(std::min(p1(0),p2(0)),std::min(p1(1),p2(1)),std::min(p1(2),p2(2)));
    this->max = Eigen::Vector3f(std::max(p1(0),p2(0)),std::max(p1(1),p2(1)),std::max(p1(2),p2(2)));
  }
  bbox(Eigen::Vector3f min,float width,float height,float depth){
    this->min=min;
    this->max=min+Eigen::Vector3f(width,height,depth);
  }

  /** update bbox to include a point */
  void expand(Eigen::Vector3f p){
    if(p(0)<this->min(0)){
      this->min(0) = p(0);
    }else if(p(0)>this->max(0)){
      this->max(0) = p(0);
    }

    if(p(1)<this->min(1)){
      this->min(1) = p(1);
    }else if(p(1)>this->max(1)){
      this->max(1) = p(1);
    }

    if(p(2)<this->min(2)){
      this->min(2) = p(2);
    }else if(p(2)>this->max(2)){
      this->max(2) = p(2);
    }
  }

  void print(){
    printf("bbox is (%.3f,%.3f,%.3f)<->(%.3f,%.3f,%.3f)\n",this->min(0),this->min(1),this->min(2),this->max(0),this->max(1),this->max(2));
  }
};

/* Used when computing tangent planes on the surface */
struct Plane{
  Eigen::Vector3f center, normal;
};

/* Defined as connection between two vertices */
struct Edge{
  int v1,v2;
  float weight;
  Edge(){
    this->v1=0;
    this->v2=0;
    this->weight=0;
  }
  Edge(int v1,int v2){
    this->v1=v1;
    this->v2=v2;
    this->weight=0;
  }
  Edge(int v1,int v2,float w){
    this->v1=v1;
    this->v2=v2;
    this->weight=w;
  }
  bool operator< (const Edge e) const{
   return (this->weight)<(e.weight);
 }
};

/* Used to group data in point clouds based on proximity */
struct CubeData{
  bbox cube;
  std::vector<Eigen::Vector3f>* vertices;
  std::vector<int>* indices;
  CubeData(Eigen::Vector3f min,Eigen::Vector3f max){
    cube = bbox(min,max);
    vertices = new std::vector<Eigen::Vector3f>(0);
    indices = new std::vector<int>(0);
  }
};

inline void printPoint(Eigen::Vector3f p){
  printf("(%.3f,%.3f,%.3f)\n",p(0),p(1),p(2));
}
