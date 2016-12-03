/* needed for vectors */
#include <vector>

/* needed for malloc and printf and scanf and what not */
#include <stdio.h>
#include <stdlib.h>

/* sqrt and square */
#include <cmath>

/* data structures */
#include "structures.h"

/* linear algebra library */
#include <util/eigen-eigen-26667be4f70b/Eigen/Dense>

using namespace std;
using namespace Eigen;

Plane getTangentPlane(vector<Vector3f> neighbors){
  Plane tangentPlane;
  tangentPlane.center(0.f,0.f,0.f);
  tangentPlane.normal(0.f,0.f,1.f);
  int numNeighbors = neighbors.size();

  MatrixXf matPoints(numNeighbors,3);

  // get centroids of neigbors
  printf("%d neighbors\n",numNeighbors);
  for(int i=0;i<numNeighbors;i++){
    tangentPlane.center += neighbors[i];
    matPoints.row(i) = neighbors[i];
  }

  tangentPlane.center/= numNeighbors;

  // covariance matrix
  matPoints -= tangentPlane.center;
  Matrix3f cov;
  cov = matPoints * matPoints.transpose();
  JacobiSVD<Matrix3f> svd(cov, ComputeThinV);

  // S is stored in svd.singularValues()
  // By default S is sorted in decreasing order. We need the smallest eigenvalue.
  // V is stored in svd.matrixV()
  // Set the Normal to the Eigenvector corresponding to the smallest eigenvalue.
  tangentPlane.normal = svd.matrixV.col(2);
  
  return tangentPlane;
}

Plane* computeTangentPlanes(Vector3f* points, int numPoints, float ro, float delta){
  Plane* planes = (Plane *) malloc(sizeof(Plane)*numPoints);
  for (int i=0;i<numPoints;i++){
    vector<Vector3f> neighbors;
    Vector3f p1 = points[i];
    for (int j=0;j<numPoints;j++){
      // if point within ro+delta, add to nearby_points
      Vector3f p2 = points[j];
      float dist = (p2-p1).norm();
      if (dist<=(ro+delta)) neighbors.push_back(points[j]);
    }
    planes[i] = getTangentPlane(neighbors);
  }
  return planes;
}

float getDist(Vector3f p, Plane* planes, int numPlanes){
  //given a point, approximate its distance to the nearest plane
  
  //first find closest plane
  int closestIdx = 0;
  float closestDist = (p-planes[0].center).norm();
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
