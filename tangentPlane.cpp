/* needed for vectors */
#include <vector>

/* needed for malloc and printf and scanf and what not */
#include <stdio.h>
#include <stdlib.h>

/* sqrt and square */
#include <cmath>

/* data structures */
#include <structures.h>

/* armadillo is a linear algebra library */
#include <armadillo>


/* computeTangentPlanes */
std::vector<TPlane> computeTangentPlanes(vector<Point>* points, float ro, float delta){
  // For each point, group points within ro+delta into nearby_points 
  std::vector<TPlane> tangent_planes (*points.size);
  int i,j,k;
  for (i=0;i<numPoints;i++){
    std::vector<Point> nearby_points;
    k=0;
    for (j=0;j<numPoints;j++){
      // if point within ro+delta, add to nearby_points
      if (dist(points[i],points[j])<=(ro+delta)){
        nearby_points[k] = points[j];
        k++;
      }
    }
    tangent_planes[i] = leastSqTPlane(nearby_points,k);
  }
}

/* leastSqTPlane */
TPlane leastSqTPlane(vector<Point>* nearby_points, int numPoints){
  TPlane tangentPlane;
  mat matrixPoints = zeros<mat>(3,numPoints);
  mat centroid = zeros<mat>(3,numPoints);
  int i;
  for(i=0;i<nearby_points.size;i++){
    // centroid computation
    tangetPlane.center.x += nearby_points[i].x;
    tangetPlane.center.y += nearby_points[i].y;
    tangetPlane.center.z += nearby_points[i].z;
    
    centroid(i,0) += nearby_points[i].x;
    centroid(i,1) += nearby_points[i].y;
    centroid(i,2) += nearby_points[i].z;

    // place points into matrix for normal computation
    matrixPoints(i,0) = nearby_points[i].x;
    matrixPoints(i,1) = nearby_points[i].y;
    matrixPoints(i,2) = nearby_points[i].z;
  }
  // centroid computation
  tangentPlane.center.x /= numPoints;
  tangentPlane.center.y /= numPoints;
  tangentPlane.center.z /= numPoints;
  centroid /= numPoints;

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
  
  TPlane.normal.x = eigvec(0,i);
  TPlane.normal.y = eigvec(1,i);
  TPlane.normal.z = eigvec(2,i);

}

