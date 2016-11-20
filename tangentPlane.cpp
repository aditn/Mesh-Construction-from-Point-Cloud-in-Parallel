/* needed for vectors */
#include <vector>

/* needed for malloc and printf and scanf and what not */
#include <stdio.h>
#include <stdlib.h>

/* sqrt and square */
#include <cmath>

/* data structures */
#include "structures.h"

/* armadillo is a linear algebra library */
//#include <armadillo>

/* svd3.h*/
#include "util/svd3/svd3.h"

/* matrixFunctions */
#include "util/matrixFunc.h"

using namespace std;

// getTangentPlane returns a tangent plane of a single point
Plane getTangentPlane(vector<V3> neighbors){
  Plane tangentPlane;
  tangentPlane.center = V3(0.f,0.f,0.f);
  tangentPlane.normal = V3(0.f,0.f,0.f);
  
  int numNeighbors = neighbors.size();
  printf("%d neighbors\n",numNeighbors);

  // init matrixNeighbors matrix
  double neighMat[numNeighbors][3];
  mat matrixNeighbors = mat(numNeighbors, 3, &neighMat);

  for(int i=0;i<numNeighbors;i++){
    tangentPlane.center.add(neighbors[i]);

    // place points into matrix for normal computation
    matrixNeighbors[i][0] = neighbors[i].x;
    matrixNeighbors[i][1] = neighbors[i].y;
    matrixNeighbors[i][2] = neighbors[i].z;
  }
  
  // centroid computation
  tangentPlane.center.scale(1.0/numNeighbors);

  // covariance matrix
  // TODO: write subtractCentroid
  subtractCentroid(matrixNeighbors, tangentPlane.center);
  
  // TODO: write matTranspose
  double neighMatT[3][numNeighbors];
  mat matrixNeighborsT = mat(3, numNeighbors, &neighMatT);
  transpose(matrixNeighbors, matrixNeighborsT);

  // Matrix Multiplication
  mat cov;
  cov.row = 3;
  cov.col = 3;
  cov.matrix[cov.row][cov.col];
  matMult(matrixPoints, matrixPointsT, cov);

  double[3][3] U,S,V;
  // get svd
  svd3((double *)U, (double *)S, (double *)V, (double *)cov);

  // eigenvector with smallest eigenvalue is the normal
  int smallEigValInd;
  double minEigVal = min(V[0][0], V[1][1]);
  minEigVal = min(minEigVal, V[2][2]);

  if (minEigVal == V[0][0]) smallEigValInd = 0;
  else if (minEigVal == V[1][1]) smallEigValInd = 1;
  else smallEigValInd = 2;
  
  tangentPlane.normal = V3(eigvec[i,smallEigValInd], eigvec[i,smallEigValInd], eigvec[i,smallEigValInd]);
  return tangentPlane;
}

// computeTangetPlanes gets the tangent planes at each point
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


/**************************************
 * Consistent Tangent Plane Orientation
 **************************************/

Plane* getNeighborPlanesWithCost(Plane* allTangentPlanes, int currentPlaneIndex, float ro, float delta){
  int numPlanes;
  Plane currentPlane = allTangentPlanes[currentPlaneIndex];
  V3 currentCentroid = currentPlane.center;
  std::vector<PlaneGraph> neighborPlanes;

  int j = 0;
  for (int i=0; i<numPlanes; i++){
    Plane possibleNeighbor = allTangentPlanes[i];
    if (currentCentroid.dist(possibleNeighbor.center) <= (delta+ro)){
      // add to Neighbor Planes if centroids are too close
      neighborPlanes.push_back(possibleNeighbor);
    }
  }

  return neighborPlanes;
}