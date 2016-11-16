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
#include <armadillo>

using namespace std;

Plane getTangentPlane(vector<V3> neighbors){
  Plane tangentPlane;
  tangentPlane.center = V3(0.f,0.f,0.f);
  tangentPlane.normal = V3(0.f,0.f,1.f);
  int numNeighbors = neighbors.size();
  printf("%d neighbors\n",numNeighbors);
  for(int i=0;i<numNeighbors;i++){
    tangentPlane.center.add(neighbors[i]);
    
    centroid(i,0) += neighbors[i].x;
    centroid(i,1) += neighbors[i].y;
    centroid(i,2) += neighbors[i].z;
    
    // place points into matrix for normal computation
    matrixPoints(i,0) = neighbors[i].x;
    matrixPoints(i,1) = neighbors[i].y;
    matrixPoints(i,2) = neighbors[i].z;
  }
  
  // centroid computation
  tangentPlane.center.scale(1.0/numNeighbors);
  centroid /= numNeighbors;

  // covariance matrix
  matrixPoints -= centroid;
  mat matrixPointsT = matrixPoints.t();
  mat cov = matrixPoints * matrixPointsT;

  fvec eigval;
  fmat eigvec;
  eig_sym(eigval, eigvec, cov);

  // eigenvector with smallest eigenvalue is the normal
  float smallEigVal=100;
  int smallEigValInd;
  for (i =0; i<3; i++){
    if (eigval[i] < smallEigVal){
      smallEigVal = eigval[i];
      smallEigValInd = i;
    }
  }
  
  tangentPlane.normal = V3(eigvec(0,i), eigvec(1,i), eigvec(2,i));
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

