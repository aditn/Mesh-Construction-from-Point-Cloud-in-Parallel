/* needed for vectors */
#include <vector>

/* needed for malloc and printf and scanf and what not */
#include <stdio.h>
#include <stdlib.h>

/* sqrt and square */
#include <cmath>

/* data structures */
#include "structures.h"



/**************************************
 * Consistent Tangent Plane Orientation
 **************************************/

Plane* getNeighborPlanes(Plane* allTangentPlanes, int currentPlaneIndex, float ro, float delta){
  int numPlanes;
  Plane currentPlane = allTangentPlanes[currentPlaneIndex];
  V3 currentCentroid = currentPlane.center;
  std::vector<Plane> neighborPlanes;

  for (int i=0; i<numPlanes; i++){
    Plane possibleNeighbor = allTangentPlanes[i];
    if (currentCentroid.dist(possibleNeighbor.center) <= (delta+ro)){
      // add to neighborPlanes if centroids are close
      neighborPlanes.push_back(possibleNeighbor);
    }
  }
  return neighborPlanes;
}

// Take neighbor planes and construct graph via kruskal algo.
// Cost function based on dot product of normals
PlaneGraph* constructPlaneGraph(Plane* neighborPlanes){
  std::vector<PlaneGraph> graphOfNeighborPlanes;
  return graphofNeighborPlanes;
}

// propogate normals


// main function call
Plane* consistentTangentPlaneOrientation(Plane* allTangentPlanes){
  for (int i=0; i<numPlanes; i++){
    
  }
}
