/* needed for vectors */
#include <vector>

/* needed for malloc and printf and scanf and what not */
#include <stdio.h>
#include <stdlib.h>

/* sqrt and square */
#include <cmath>

/* data structures */
#include <structures.h>

/* computeTangentPlanes */
std::vector<TPlane> computeTangentPlanes(vector<Point>* points, float ro, float delta){
  // For each point, group points within ro+delta into nearby_points 
  std::vector<TPlane> tangent_planes (*points.size);
  int i,j,k;
  for (i=0;i<numPoints;i++){
    std::vector<Point> nearby_points;
    k=0;
    for (j=0;j<numPoints;j++){
      // if point within ro+delta, add to nearby_points
      if (dist(points[i],points[j])<=(ro+delta)){
        nearby_points[k] = points[j];
        k++;
      }
    }
    tangent_planes[i] = leastSqTPlane(nearby_points,k);
  }
}

/* leastSqTPlane */
TPlane leastSqTPlane(vector<Point>* nearby_points, int numPoints){
  TPlane tangentPlane;
  int i;
  for(i=0;i<nearby_points.size;i++){
    tangetPlane.center.x += nearby_points[i].x;
    tangetPlane.center.y += nearby_points[i].y;
    tangetPlane.center.z += nearby_points[i].z;
  }
  tangentPlane.center.x /= numPoints;
  tangentPlane.center.y /= numPoints;
  tangentPlane.center.z /= numPoints;
}
