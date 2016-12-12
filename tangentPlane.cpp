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

Plane* computeTangentPlanes(Vector3f* points, int numPoints){
  Plane* planes = (Plane *) malloc(sizeof(Plane)*numPoints);
  
  // Get bounding size of entire point cloud
  Vector3f minxyz(FLT_MAX,FLT_MAX,FLT_MAX),
           maxxyz(FLT_MIN,FLT_MIN,FLT_MIN);
  for (int i=0; i<numPoints; i++){
    if (points[i](0)<minxyz(0)) minxyz(0) = points[i](0);
    if (points[i](1)<minxyz(1)) minxyz(1) = points[i](1);
    if (points[i](2)<minxyz(2)) minxyz(2) = points[i](2);
    if (points[i](0)>maxxyz(0)) maxxyz(0) = points[i](0);
    if (points[i](1)>maxxyz(1)) maxxyz(1) = points[i](1);
    if (points[i](2)>maxxyz(2)) maxxyz(2) = points[i](2);
  }
  Vector3f diff = maxxyz - minxyz;
  
  // increase of universe to allow for even number of cubes
  float cubeLength = (rho+delta);
  float incrX = int(ceil(diff(0)/cubeLength)+1)*cubeLength - diff(0);
  float incrY = int(ceil(diff(1)/cubeLength)+1)*cubeLength - diff(1);
  float incrZ = int(ceil(diff(2)/cubeLength)+1)*cubeLength - diff(2);
  Vector3f increase(incrX,incrY,incrZ);
  maxxyz += increase/2; 
  minxyz -= increase/2;
  diff = maxxyz - minxyz;
  printPoint(diff);
 
  Vector3f numCubesDim = diff/cubeLength;
  printPoint(numCubesDim);
  int numCubes = int(numCubesDim(0)*numCubesDim(1)*numCubesDim(2));
  printf("numCubes:%d, %f\n",numCubes, numCubesDim(0)*numCubesDim(1)*numCubesDim(2));
  // Define Array of CubeData
  CubeData*** splitData = (CubeData ***) malloc(numCubesDim(0)*sizeof(CubeData**));
  for (int i=0; i<int(numCubesDim(0)); i++){
    splitData[i] = (CubeData **) malloc(numCubesDim(1)*sizeof(CubeData*));
    for (int j=0; j<int(numCubesDim(1)); j++){
      splitData[i][j] = (CubeData *) malloc(numCubesDim(2)*sizeof(CubeData));
      for (int k=0; k<int(numCubesDim(2)); k++){
        CubeData cd;
        cd.cube = bbox(minxyz+Vector3f(i*cubeLength,j*cubeLength,k*cubeLength),
                       minxyz+Vector3f((i+1)*cubeLength,(j+1)*cubeLength,(k+1)*cubeLength));
        //vector<Vector3f> v;
        //cd.vertices = v;
        splitData[i][j][k] = cd;
        //printf("i:%d,j:%d,k:%d\n", i,j,k);
      }
    }
  } 



  /*CubeData* splitData = (CubeData *) calloc(int(numCubes), sizeof(CubeData));
  //Vector3f curCoord = minxyz;
  //Vector3f prevCoord = curCoord;
  //Vector3f incr(cubeLength,cubeLength,cubeLength);
  printf("made space for cubes\n");
  // Set bbox for each Cube
  int f = 0;
  for(float z=minxyz(2); z<maxxyz(2); z+=cubeLength){
    for(float y=minxyz(1); y<maxxyz(1); y+=cubeLength){
      for(float x=minxyz(0); x<maxxyz(0); x+=cubeLength){
        splitData[f].cube = bbox(Vector3f(x,y,z),Vector3f(x+cubeLength,y+cubeLength,z+cubeLength));
        f+=1;
      }
    }
  }*/

  //printf("a: \n");
  //printPoint(splitData[106].cube.min);
  //printPoint(splitData[106].cube.max);
  printf("assigned bbox info for each cube\n");
  printf("maxxyz\n");
  printPoint(maxxyz);
  // Put points into specific CubeData struct (O(n))
  insertPoints(points, numPoints, splitData, maxxyz, minxyz, diff, cubeLength, numCubes);

  printf("placed points into cube struct\n");
  

  //int width = int(diff(1)/cubeLength);
  //int depth = int(diff(2)/cubeLength);
  
  // Compute tangent planes by looking within specific cubes
  for (int i=0; i<numPoints; i++){
    vector<Vector3f> neighbors;
    Vector3f curPoint = points[i];
    Vector3f cubesToSearchInd[27];
    //Vector3f cubecoord = (curPoint-minxyz)/cubeLength;

    int indA=0;
    int indB=0;
    int indC=0;
    Vector3f cubeInd;
    for (int a=-1; a<2; a++){ // determine which cubes to search through
      for (int b=-1; b<2; b++){
        for (int c=-1; c<2; c++){
          //printf("indC: %d\n", indC);
          Vector3f neighPoint(curPoint(0)+a*cubeLength,
                              curPoint(1)+b*cubeLength,
                              curPoint(2)+c*cubeLength);
          if ((neighPoint(0)>maxxyz(0) || neighPoint(0)<minxyz(0)) ||
              (neighPoint(1)>maxxyz(1) || neighPoint(1)<minxyz(1)) ||
              (neighPoint(2)>maxxyz(2) || neighPoint(2)<minxyz(2))){
            cubesToSearchInd[indC] = Vector3f(-1,-1,-1);
            indC += 1;
            continue;
          }
          cubeInd = (neighPoint-minxyz)/cubeLength;
          /*cubesToSearchInd[ind] = int(cubeInd(0))*width*depth+
                                  int(cubeInd(1))*depth+
                                  int(cubeInd(2));*/
          cubesToSearchInd[indC] = Vector3f(int(cubeInd(0)),int(cubeInd(1)),int(cubeInd(2)));
          indC += 1;
        }
        indB += 3;
        indC = indB;
      }
      indA += 9;
      indB = indA;
    }

    //printf("cubeInd:\n");
    //printf("got cubes to look through\n");

    for(int k=0; k<27; k++){ // go through list of each cube
      Vector3f cubeIndex = cubesToSearchInd[k];
      //printf("cubeIndex:");
      //printPoint(cubeIndex);
      if (cubeIndex(0) == -1) continue;
      CubeData c = splitData[int(cubeIndex(0))][int(cubeIndex(1))][int(cubeIndex(2))];
      if (c.vertices.size()==0) continue;
      for (unsigned j=0; j<c.vertices.size(); j++){
        //printf("cubeIndex: ");
        //printPoint(cubeIndex);
        Vector3f p2 = c.vertices[j];
        float dist = (p2-curPoint).norm();
        if (dist<=(rho+delta)) neighbors.push_back(p2);
      }
    }
    planes[i] = getTangentPlane(neighbors);
  }
  /*#ifdef USE_OMP
  int cubeIndex;
  #pragma omp parallel for
  for (int i=0; i<numPoints; i++){
    vector<Vector3f> neighbors;
    Vector3f curPoint = points[i]/cubeLength; // normalize to fit cube dim
    cubeIndex = curPoint(0)+diff(0)*(curPoint(1)+diff(2)*curPoint(2)); // maps 3d to 1d
    for (unsigned j=0; j<splitData[cubeIndex].vertices.size(); j++){
      Vector3f p2 = splitData[cubeIndex].vertices[j];
      float dist = (p2-curPoint).norm();
      if (dist<=(rho+delta)) neighbors.push_back(p2);
    }
    planes[i] = getTangentPlane(neighbors);
  }
  
  #else
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
  #endif
  */

  /*for (int i=0; i<numPoints; i++){
    printf("tangent%d, centroid, normal:",i);
    printPoint(planes[i].center);
    printPoint(planes[i].normal);
  }*/


  return planes;
}

/*#ifdef USE_OMP
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
#endif*/




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

void insertPoints(Eigen::Vector3f* points, int numPoints, CubeData*** splitData,
  Eigen::Vector3f maxxyz,Eigen::Vector3f minxyz, Eigen::Vector3f diff, float cubeLength, int numCubes){
  
  Vector3f curPoint, cubecoord, curInd;
  //int cubeInd;
  //int width = int(diff(1)/cubeLength);
  //int depth = int(diff(2)/cubeLength);
  for(int i=0; i<numPoints; i++){
    curPoint = points[i];
    cubecoord = (curPoint-minxyz)/cubeLength;
    //printPoint(cubecoord);
    /*cubeInd = int(cubecoord(0))*width*depth+
              int(cubecoord(1))*depth+
              int(cubecoord(2));*/
    curInd = Vector3f(floor(cubecoord(0)),floor(cubecoord(1)),floor(cubecoord(2)));
    //cubeInd = int(cubecoord(2)) + width*(int(cubecoord(1))+depth*int(cubecoord(0)));
    /*printf("quantPoint: ");
    printPoint(cubecoord);
    printf("curInd: ");
    printPoint(curInd);
    */
    Vector3f bboxmin = splitData[int(curInd(0))][int(curInd(1))][int(curInd(2))].cube.min;
    Vector3f bboxmax = splitData[int(curInd(0))][int(curInd(1))][int(curInd(2))].cube.max;
    
    /*printf("cubeInd: ");
    printPoint(curInd);
    printf("bbox dims (min, max)\n" );
    printPoint(bboxmin);
    printPoint(bboxmax);
    printf("pointCur: ");
    printPoint(curPoint);
    */
    if ((curPoint(0)<bboxmin(0) || curPoint(0)> bboxmax(0)) ||
        (curPoint(1)<bboxmin(1) || curPoint(1)> bboxmax(1)) ||
        curPoint(2)<bboxmin(2) || curPoint(2)> bboxmax(2)){
      // point shouldn't be in this cube
      printf("cubeInd: ");
      printPoint(curInd);
      printf("bbox dims (min, max)\n" );
      printPoint(bboxmin);
      printPoint(bboxmax);
      printf("pointCur: ");
      printPoint(curPoint);
    }
    //printf("booya\n");
    /*if (cubeInd>numCubes){
      printf("cubeIndold: %d\n",cubeInd);
      printf("coordOld: ");
      printPoint(cubecoord+minxyz);
      for (int i = 0; i < 3; i++){
        if (cubecoord(i)+minxyz(i) >= maxxyz(i)){
          cubecoord(i)-=cubeLength;
        }
        else if (cubecoord(i)+minxyz(i) <= minxyz(i)){
          cubecoord(i)+=cubeLength;
        }
      }
      cubeInd = int(cubecoord(2)) + width*(int(cubecoord(1))+depth*int(cubecoord(0)));
      printf("cubeIndnew: %d\n",cubeInd);
      printf("coordNew: ");
      printPoint(cubecoord+minxyz);
    }
    if (cubeInd==106) {
      printf("hello");
      printPoint(points[i]);
    }*/
    //printf("curInd:%d, %d, %d\n",int(curInd(0)),int(curInd(1)),int(curInd(2)));
    CubeData c = splitData[int(curInd(0))][int(curInd(1))][int(curInd(2))];
    //printf("can access cube\n");
    c.vertices.push_back(curPoint);
  }
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
