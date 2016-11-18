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

  //TODO: compute tangent
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
