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
  void print(){
    printf("(%.3f,%.3f,%.3f)\n",this->x,this->y,this->z);
  }
};
inline V3 operator-(const V3& a, const V3& b){
    V3 out;
    out.x = a.x-b.x;
    out.y = a.y-b.y;
    out.z = a.z-b.z;
    return out;
}
inline V3 operator+(const V3& a, const V3& b){
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
inline V3 operator*(const V3& a,const float& b){
  V3 out = V3(a);
  out.scale(b);
  return out;
}
inline V3 operator*(const float& b,const V3& a){
  V3 out = V3(a);
  out.scale(b);
  return out;
}

struct bbox{
  V3 min,max;
  bbox(){
    this->min = V3();
    this->max = V3();
  }
  bbox(V3 p1, V3 p2){
    this->min = V3(std::min(p1.x,p2.x),std::min(p1.y,p2.y),std::min(p1.z,p2.z));
    this->max = V3(std::max(p1.x,p2.x),std::max(p1.y,p2.y),std::max(p1.z,p2.z));
  }
  bbox(V3 min,float width,float height,float depth){
    this->min=V3(min);
    this->max=min+V3(width,height,depth);
  }

  /** update bbox to include a point */
  void expand(V3 p){
    if(p.x<this->min.x){
      this->min.x = p.x;
    }else if(p.x>this->max.x){
      this->max.x = p.x;
    }

    if(p.y<this->min.y){
      this->min.y = p.y;
    }else if(p.y>this->max.y){
      this->max.y = p.y;
    }

    if(p.z<this->min.z){
      this->min.z = p.z;
    }else if(p.z>this->max.z){
      this->max.z = p.z;
    }
  }

  void print(){
    printf("bbox is (%.3f,%.3f,%.3f)<->(%.3f,%.3f,%.3f)\n",this->min.x,this->min.y,this->min.z,this->max.x,this->max.y,this->max.z);
  }
};

struct Plane{
  V3 center, normal;
};

struct Edge{
  V3 v1,v2;
  Edge(){
  this->v1 = V3();
  this->v2 = V3();
  }
  Edge(V3 v1,V3 v2){
    this->v1 = V3(v1);
    this->v2 = V3(v2);
  }
};

struct E{
  int v1,v2;
  E(){
    this->v1=0;
    this->v2=0;
  }
  E(int v1,int v2){
    this->v1=v1;
    this->v2=v2;
  }
};

inline void printPoint(V3 p){
  printf("(%.3f,%.3f,%.3f)\n",p.x,p.y,p.z);
}
