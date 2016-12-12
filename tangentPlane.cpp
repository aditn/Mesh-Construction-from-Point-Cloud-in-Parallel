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
    if(neighbors.size()<3) printf("neighbors for %d are %d\n", i,(int)neighbors.size());
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
      if (dist<=(rho+delta)) neighbors.push_back(points[j]);
    }
    planes[i] = getTangentPlane(neighbors);
  }
  return planes;
}
#endif




inline void insMin(int* idxs, float* vals, int N, int newIdx, float newVal){
  int curIdx = N-1;

  if(vals[curIdx]>newVal){
    vals[curIdx] = newVal;
    idxs[curIdx] = newIdx;
    curIdx--;
  }else return;

  while(curIdx>=0 && vals[curIdx]>newVal){
    vals[curIdx+1] = vals[curIdx];
    idxs[curIdx+1] = idxs[curIdx];
    vals[curIdx] = newVal;
    idxs[curIdx] = newIdx;
    curIdx--;
  }
}

std::vector<int> getNearest(Vector3f* points, int numPoints, int idx, int numNeighbors){
  std::vector<int> neighbors;
  if(numPoints<=numNeighbors){
    for(int i=0;i<numPoints;i++) neighbors.push_back(i);
    return neighbors;
  }else{
    int* nearestIdxs = (int*) malloc(sizeof(int)*numNeighbors);
    float* nearestVals = (float*) malloc(sizeof(float)*numNeighbors);
    for(int i=0;i<numNeighbors;i++) nearestVals[i] = INFINITY;

    Vector3f curPoint = points[idx];
    for(int i=0;i<numPoints;i++){
      if(i==idx) continue;
      insMin(nearestIdxs,nearestVals,numNeighbors,i,(curPoint-points[i]).norm());
    }
    for(int i=0;i<numNeighbors;i++) neighbors.push_back(nearestIdxs[i]);
    free(nearestIdxs);
    free(nearestVals);
    return neighbors;
  }
}







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

void setDists(float* dists, Vector3f* ps, int numPoints, Plane* planes, int numPlanes){
  int* closestIdxs = (int*)   calloc(numPoints,sizeof(int));

  for(int j=0;j<numPoints;j++){ //assume plane0 is the closest
    dists[j] = (ps[j]-planes[0].center).norm();
  }

  for(int i=1;i<numPlanes;i++){
    Vector3f cur = planes[i].center;
    for(int j=0;j<numPoints;j++){
      float curDist = (ps[j]-cur).norm();
      if(curDist<dists[j]){
        dists[j]=curDist;
        closestIdxs[j] = i;
      }
    }
  }

  for(int i=0;i<numPoints;i++){
    Plane curPlane = planes[closestIdxs[i]];
    dists[i] = (ps[i]-curPlane.center).dot(curPlane.normal);
  }
}
