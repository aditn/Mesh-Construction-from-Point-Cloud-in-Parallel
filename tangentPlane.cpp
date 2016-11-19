/* needed for vectors */
#include <vector>

/* needed for malloc and printf and scanf and what not */
#include <stdio.h>
#include <stdlib.h>

/* sqrt and square */
#include <cmath>

/* data structures */
#include "structures.h"

using namespace std;

Plane getTangentPlane(vector<V3> neighbors){
  Plane tangentPlane;
  tangentPlane.center = V3(0.f,0.f,0.f);
  tangentPlane.normal = V3(0.f,0.f,1.f);
  int numNeighbors = neighbors.size();
  for(int i=0;i<numNeighbors;i++) tangentPlane.center.add(neighbors[i]);
  tangentPlane.center.scale(1.0/numNeighbors);

  //TODO: compute tangent the right way
  if(numNeighbors>=3){
    V3 v1 = neighbors[1]-neighbors[0];
    V3 v2 = neighbors[2]-neighbors[0];
    tangentPlane.normal = v1.cross(v2);
  }

  return tangentPlane;
}

Plane* computeTangentPlanes(V3* points, int numPoints, float ro, float delta){
  Plane* planes = (Plane *) malloc(sizeof(Plane)*numPoints);
  for (int i=0;i<numPoints;i++){
    vector<V3> neighbors;
    V3 p1 = points[i];
    for (int j=0;j<numPoints;j++){
      // if point within ro+delta, add to nearby_points
      V3 p2 = points[j];
      float dist = p1.dist(p2);
      if (dist<=(ro+delta)) neighbors.push_back(points[j]);
    }
    planes[i] = getTangentPlane(neighbors);
  }
  return planes;
}

float getDist(V3 p, Plane* planes, int numPlanes){
  //given a point, approximate its distance to the nearest plane
  
  //first find closest plane
  int closestIdx = 0;
  float closestDist = p.dist(planes[0].center);
  for(int i=1;i<numPlanes;i++){
    float curDist = p.dist(planes[i].center);
    if(curDist<closestDist){
      closestDist = curDist;
      closestIdx = i;
    }
  }

  //then approximate dist to surface
  return (p-planes[closestIdx].center).dot(planes[closestIdx].normal); 
}
