/* needed for vectors */
#include <vector>

/* needed for malloc and printf and scanf and what not */
#include <stdio.h>
#include <stdlib.h>

/* sqrt and square */
#include <cmath>

/* data structures */
#include "structures.h"
#include "tangentPlane.h"
#include "constants.h"

/* linear algebra library */
#include <Eigen/Dense>

#include <omp.h>

using namespace std;
using namespace Eigen;

Plane getTangentPlane(vector<Vector3f> neighbors){
  Plane tangentPlane;
  tangentPlane.center = Vector3f(0.f,0.f,0.f);
  tangentPlane.normal = Vector3f(0.f,0.f,0.f);
  int numNeighbors = neighbors.size();

  MatrixXf matPoints(numNeighbors,3);

  // get centroids of neigbors
  for(int i=0;i<numNeighbors;i++){
    tangentPlane.center += neighbors[i];
    matPoints.row(i) = neighbors[i];
  }

  tangentPlane.center/= numNeighbors;

  // covariance matrix
  for(int i=0;i<numNeighbors;i++){
    matPoints.row(i) -= tangentPlane.center;
  }

  MatrixXf cov = MatrixXf::Zero(3,3);
  cov = matPoints.transpose() * matPoints;
  JacobiSVD<MatrixXf> svd(cov, ComputeThinV);

  // S is stored in svd.singularValues()
  // By default S is sorted in decreasing order. We need the smallest eigenvalue.
  // V is stored in svd.matrixV()
  // Set the Normal to the Eigenvector corresponding to the smallest eigenvalue.
  tangentPlane.normal = svd.matrixV().col(2);
  
  return tangentPlane;
}

#ifdef USE_OMP
Plane* computeTangentPlanes(Vector3f* points, int numPoints){
  Plane* planes = (Plane *) malloc(sizeof(Plane)*numPoints);

  bool* close_enough = (bool*) calloc(numPoints*numPoints,sizeof(bool));

  #pragma omp parallel for
  for(int i=0;i<numPoints*numPoints;i++){
    int row = i/numPoints,
        col = i%numPoints;
    float dist = (points[col]-points[row]).norm();
    if(dist <= (rho+delta)) close_enough[i] = true;
  }

  #pragma omp parallel for
  for(int i=0;i<numPoints;i++){
    vector<Vector3f> neighbors;
    int idx_coord = numPoints*i;
    for(int idx=0;idx<numPoints;idx++){
      if(close_enough[idx_coord+idx]) neighbors.push_back(points[idx]);
    }
    planes[i] = getTangentPlane(neighbors);
  }

  return planes;
}
#else
Plane* computeTangentPlanes(Vector3f* points, int numPoints){
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
#endif

float getDist(Vector3f p, Plane* planes, int numPlanes){
  //given a point, approximate its distance to the nearest plane
  
  //first find closest plane
  int closestIdx = 0;
  float closestDist = (p-planes[0].center).norm();
  for(int i=1;i<numPlanes;i++){
    float curDist = (p-planes[i].center).norm();
    if(curDist<closestDist){
      closestDist = curDist;
      closestIdx = i;
    }
  }
  //then approximate dist to surface
  return (p-planes[closestIdx].center).dot(planes[closestIdx].normal); 
}
