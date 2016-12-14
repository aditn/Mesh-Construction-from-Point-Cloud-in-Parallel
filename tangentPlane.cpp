/* needed for vectors */
#include <vector>

/* needed for malloc and printf and scanf and what not */
#include <stdio.h>
#include <stdlib.h>

/* sqrt and square */
#include <cmath>
#include <float.h>

/* data structures */
#include "structures.h"
#include "tangentPlane.h"
#include "constants.h"
#include "ourTimer.h"
#include "tangentPlane.h"

/* linear algebra library */
#include <Eigen/Dense>

#include <omp.h>

using namespace std;
using namespace Eigen;

extern float rho;
CubeData*** LookupStructure;
int dimx,dimy,dimz;
bbox universe;
float cubeLength;

Plane getTangentPlane(vector<Vector3f> neighbors){
  Plane tangentPlane;
  tangentPlane.center = Vector3f(0.f,0.f,0.f);
  tangentPlane.normal = Vector3f(0.f,0.f,0.f);
  int numNeighbors = neighbors.size();

  MatrixXf matPoints(numNeighbors,3);

  // get centroids of neighbors
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

void insertPoints(Eigen::Vector3f* points, int numPoints, CubeData*** splitData,bbox sysBox,float cubeLength){
  for(int i=0; i<numPoints; i++){
    Vector3f curPoint = points[i];
    Vector3f cubecoord = (curPoint-sysBox.min)/cubeLength;
    int x=floor(cubecoord(0)),
        y=floor(cubecoord(1)),
        z=floor(cubecoord(2));
    splitData[x][y][z].vertices->push_back(curPoint);
    splitData[x][y][z].indices->push_back(i);
  }
}

Plane* computeTangentPlanes(Vector3f* points, int numPoints){
  Plane* planes = (Plane *) malloc(sizeof(Plane)*numPoints);
 
  printf("starting compute tangent Planes\n");
  printf("time: %.4fs\n",timeSince());
 
  // Get bounding size of entire point cloud
  universe=bbox();
  for(int i=0;i<numPoints;i++) universe.expand(points[i]);
  Vector3f diff = universe.max-universe.min;
  cubeLength=rho;
  // increase of universe to allow for even number of cubes
  float incrX = int(ceil(diff(0)/cubeLength)+1)*cubeLength - diff(0);
  float incrY = int(ceil(diff(1)/cubeLength)+1)*cubeLength - diff(1);
  float incrZ = int(ceil(diff(2)/cubeLength)+1)*cubeLength - diff(2);
  Vector3f increase = 0.5*Vector3f(incrX,incrY,incrZ);
  universe.max += increase; 
  universe.min -= increase;
  Vector3f numCubesDim = (universe.max-universe.min)/cubeLength;
  printf("Num cubes in each dir is:\n");
  printPoint(numCubesDim);
  dimx=int(numCubesDim(0));
  dimy=int(numCubesDim(1));
  dimz=int(numCubesDim(2));
  printf("Using %dx%dx%d\n",dimx,dimy,dimz);
  // Define Array of CubeData
  CubeData*** splitData = (CubeData ***) malloc(dimx*sizeof(CubeData**));
  for (int i=0; i<dimx; i++){
    splitData[i] = (CubeData **) malloc(dimy*sizeof(CubeData*));
    for (int j=0; j<dimy; j++){
      splitData[i][j] = (CubeData *) malloc(dimz*sizeof(CubeData));
      for (int k=0; k<dimz; k++){
        splitData[i][j][k] = CubeData(universe.min+Vector3f(i*cubeLength,j*cubeLength,k*cubeLength),universe.min+Vector3f((i+1)*cubeLength,(j+1)*cubeLength,(k+1)*cubeLength));
      }
    }
  }
  printf("structure allocated.\n");
  printf("time: %.4fs\n",timeSince());
  // Put points into specific CubeData struct (O(n))
  insertPoints(points, numPoints, splitData, universe, cubeLength);
  printf("placed points into cube structs\n");
  printf("time: %.4fs\n",timeSince());
 
  // Compute tangent planes by looking within specific cubes
#ifdef USE_OMP
  #pragma omp parallel for
#endif
  for (int pi=0; pi<numPoints; pi++){
    vector<Vector3f> neighbors;
    Vector3f curPoint = points[pi];
    int x,y,z;
    x=(curPoint(0)-universe.min(0))/cubeLength;
    y=(curPoint(1)-universe.min(1))/cubeLength;
    z=(curPoint(2)-universe.min(2))/cubeLength;

    //check nearby cubes
    for(int i=x-1;i<=x+1;i++){
      for(int j=y-1;j<=y+1;j++){
        for(int k=z-1;k<=z+1;k++){
          if(i<0 || i>=dimx ||
             j<0 || j>=dimy ||
             k<0 || k>=dimz) continue; //out of bounds
          CubeData curData = splitData[i][j][k];
          vector<Vector3f> curverts = *(curData.vertices);
          for(unsigned idx=0;idx<curverts.size();idx++){
            float dist = (curverts[idx] - curPoint).norm();
            if(dist<rho) neighbors.push_back(curverts[idx]);
          }
        }
      }
    }
    planes[pi] = getTangentPlane(neighbors);
  }

  for(int i=0;i<dimx;i++){
    for(int j=0;j<dimy;j++){
      //free vertex buf?
      free(splitData[i][j]);
    }
    free(splitData[i]);
  }

  printf("planes made!\n");
  printf("time: %.4fs\n",timeSince());
 
  return planes;
}


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

/*
vector<int> getNearest(Vector3f* points, int numPoints, int idx, int numNeighbors){
  vector<int> neighbors;
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

vector<Edge> getNeighborEdges(Vector3f* points,int numPoints,Plane* planes){
  vector<Edge> edges;
#ifdef USE_OMP
  vector<int>* neighbors = (vector<int>*) malloc(sizeof(std::vector<Edge>)*numPoints);
  #pragma omp parallel for
  for(int i=0;i<numPoints;i++){
    neighbors[i] = getNearest(points,numPoints,i,K);
  }
  for(int i=0;i<numPoints;i++){
    for(unsigned int j=0;j<neighbors[i].size();j++){
      int idx = neighbors[i][j];
      edges.push_back(Edge(i,idx,1-fabs(planes[i].normal.dot(planes[idx].normal))));
    }
  }
  free(neighbors);
#else
  for(int i=0;i<numPoints;i++){
    vector<int> neighbors = getNearest(points,numPoints,i,K);
    for(unsigned int j=0;j<neighbors.size();j++){
      int idx = neighbors[j];
      edges.push_back(Edge(i,idx,1-fabs(planes[i].normal.dot(planes[idx].normal))));
    }
  }
#endif
  return edges;
}
*/

vector<Edge> getNeighborEdges(Vector3f* points,int numPoints,Plane* planes){ 
  // Get bounding size of entire point cloud
  universe = bbox();
  universe.print();
  for(int i=0;i<numPoints;i++) universe.expand(points[i]);
  Vector3f diff = universe.max-universe.min;
  universe.print();
  cubeLength=rho;
  // increase of universe to allow for even number of cubes
  float incrX = int(ceil(diff(0)/cubeLength)+1)*cubeLength - diff(0);
  float incrY = int(ceil(diff(1)/cubeLength)+1)*cubeLength - diff(1);
  float incrZ = int(ceil(diff(2)/cubeLength)+1)*cubeLength - diff(2);
  Vector3f increase = 0.5*Vector3f(incrX,incrY,incrZ);
  universe.max += increase; 
  universe.min -= increase;
  universe.print();
  printf("side len is %.3f\n",cubeLength);
  Vector3f numCubesDim = (universe.max-universe.min)/cubeLength;
  printf("Num cubes in each dir is:\n");
  printPoint(numCubesDim);
  dimx=int(numCubesDim(0));
  dimy=int(numCubesDim(1));
  dimz=int(numCubesDim(2));
  printf("Using %dx%dx%d\n",dimx,dimy,dimz);
  // Define Array of CubeData
  LookupStructure = (CubeData ***) malloc(dimx*sizeof(CubeData**));
  for (int i=0; i<dimx; i++){
    LookupStructure[i] = (CubeData **) malloc(dimy*sizeof(CubeData*));
    for (int j=0; j<dimy; j++){
      LookupStructure[i][j] = (CubeData *) malloc(dimz*sizeof(CubeData));
      for (int k=0; k<dimz; k++){
        LookupStructure[i][j][k] = CubeData(universe.min+Vector3f(i*cubeLength,j*cubeLength,k*cubeLength),universe.min+Vector3f((i+1)*cubeLength,(j+1)*cubeLength,(k+1)*cubeLength));
      }
    }
  }
  printf("structure allocated.\n");
  printf("time: %.4fs\n",timeSince());
  // Put points into specific CubeData struct (O(n))
  insertPoints(points, numPoints, LookupStructure, universe, cubeLength);
  printf("placed points into cube structs\n");
  printf("time: %.4fs\n",timeSince());
 
  // Compute nearest neighbors by looking at nearby cubes
  vector<Edge> edges;
  vector<int>* neighbors = (vector<int>*) malloc(sizeof(std::vector<Edge>)*numPoints);
#ifdef USE_OMP
  #pragma omp parallel for
#endif
  for (int pi=0; pi<numPoints; pi++){
    vector<int> pneighbors;
    Vector3f curPoint = points[pi];
    int x,y,z;
    x=(curPoint(0)-universe.min(0))/cubeLength;
    y=(curPoint(1)-universe.min(1))/cubeLength;
    z=(curPoint(2)-universe.min(2))/cubeLength;

    int nearestIdxs[NUM_NEIGHBORS];
    float nearestVals[NUM_NEIGHBORS];
    for(int i=0;i<NUM_NEIGHBORS;i++) nearestVals[i] = INFINITY;
    int curCount = 0;

    //check nearby cubes
    bool breakNext=false;
    for(int layer=0;layer<5;layer++){ // look max 5 out
      if(curCount>NUM_NEIGHBORS){
        if(breakNext) break;
        else breakNext=true;
      }
      for(int i=x-layer;i<=x+layer;i++){
        for(int j=y-layer;j<=y+layer;j++){
          for(int k=z-layer;k<=z+layer;k++){
            if(i<0 || i>=dimx ||
               j<0 || j>=dimy ||
               k<0 || k>=dimz) continue; //out of bounds
            if(abs(i-x)<layer && abs(j-y)<layer && abs(k-z)<layer) continue; //previous layer
            CubeData curData = LookupStructure[i][j][k];
            vector<Vector3f> curverts = *(curData.vertices);
            vector<int> curindices = *(curData.indices);
            for(unsigned idx=0;idx<curverts.size();idx++){
              float dist = (curverts[idx] - curPoint).norm();
              insMin(nearestIdxs,nearestVals,NUM_NEIGHBORS,curindices[idx],dist);
              curCount++;
            }
          }
        }
      }
    }
    for(int i=0;i<min(NUM_NEIGHBORS,curCount);i++) pneighbors.push_back(nearestIdxs[i]);
    neighbors[pi] = pneighbors;
  }

  //set edges
  for(int i=0;i<numPoints;i++){
    for(unsigned int j=0;j<neighbors[i].size();j++){
      int idx = neighbors[i][j];
      edges.push_back(Edge(i,idx,1-fabs(planes[i].normal.dot(planes[idx].normal))));
    }
  }
  
  free(neighbors);

  printf("neighbor edges made!\n");
  printf("time: %.4fs\n",timeSince());
 
  return edges;
}

float getDistLookup(Vector3f p,Plane* planes,int numPlanes){ //only makes sense if querying inside bbox
  int x,y,z;
  x=(p(0)-universe.min(0))/cubeLength;
  y=(p(1)-universe.min(1))/cubeLength;
  z=(p(2)-universe.min(2))/cubeLength;
  int closestIdx = -1;
  float closestDist = 0;
  bool breakNext=false;
  int layer=0;
  while(true){ //assumes there is at least one point in the structure...
    if(closestIdx>=0){
      if(breakNext) break;
      else breakNext = true;
    }
    for(int i=x-layer;i<=x+layer;i++){
      for(int j=y-layer;j<=y+layer;j++){
        for(int k=z-layer;k<=z+layer;k++){
          if(i<0 || i>=dimx ||
             j<0 || j>=dimy ||
             k<0 || k>=dimz) continue; //out of bounds
          if(abs(i-x)<layer && abs(j-y)<layer && abs(k-z)<layer) continue; //previous layer
          CubeData curData = LookupStructure[i][j][k];
          vector<Vector3f> curverts = *(curData.vertices);
          vector<int> curindices = *(curData.indices);
          for(unsigned idx=0;idx<curverts.size();idx++){
            float dist = (curverts[idx] - p).norm();
            if(closestIdx<0 || dist<closestDist){
              closestIdx = curindices[idx];
              closestDist = dist;
            }
          }
        }
      }
    }
    layer++;
  }
  return (p-planes[closestIdx].center).dot(planes[closestIdx].normal); 
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
  for(int i=0;i<numPoints;i++) dists[i] = getDistLookup(ps[i],planes,numPlanes);
  /*
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
  }*/
}
