#define MAX_LINE_SIZE 1024 //our max line size for obj file

typedef struct{
  float x,y,z;
} Point;

typedef struct{
  Point center, normal;
} TPlane;
