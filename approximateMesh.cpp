/* needed for malloc, printf */
#include <stdio.h>
#include <stdlib.h>

/* needed for vectors */
#include <vector>
#include <algorithm> //for sorting vectors

/* fabs, square, sqrt */
#include <cmath>

/* needed for openMP! */
#include <omp.h>

/* data structures */
#include "structures.h"
#include "tangentPlane.h"
#include "parseOBJ.h"
#include "approximateMesh.h"

/* linear algebra library */
#include <Eigen/Dense>

#define K 5 //num neighbors for MST Propogation

#define USE_OMP
#ifdef USE_OMP
extern int numThreads;
#endif

extern bool DEBUG;

using namespace Eigen;

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

void approximateMesh(Vector3f* points, int numPoints, float rho, float delta,std::vector<Vector3f>& finalVertices,std::vector<Edge>& finalEdges){
  //step 1: create a plane for each point based off of its neighbors
  printf("getting planes...\n");
  printf("%d number of threads seen\n",numThreads);
  Plane* planes = computeTangentPlanes(points,numPoints,rho,delta);
  Vector3f* newPoints = (Vector3f*) malloc(sizeof(Vector3f)*numPoints);
#ifdef USE_OMP
  #pragma omp parallel for
#endif
  for(int i=0;i<numPoints;i++) newPoints[i] = planes[i].center;
  //step 2: propogate normal directions of every plane
  //substep: find plane centroid with largest z coordinate
  printf("getting max z...\n");
  int curmax = 0;
  int curidx = -1;
  for(int i=0;i<numPoints;i++){
    if(curidx==-1 || newPoints[i](2)>curmax){
      curmax = newPoints[i](2);
      curidx = i;
    }
  }
  planes[curidx].normal = Vector3f(0,0,(curmax>=0)? -1: 1); //TODO: why?

  //substep: create graph of neighbors with edge weights
  printf("getting neighbor graph...\n");
  std::vector<Edge> neighbor_edges;
  for(int i=0;i<numPoints;i++){
    std::vector<int> neighbors = getNearest(newPoints,numPoints,i,K);
    for(unsigned int j=0;j<neighbors.size();j++){
      int idx = neighbors[j];
      neighbor_edges.push_back(Edge(i,idx,1-fabs(planes[i].normal.dot(planes[idx].normal))));
    }
  }

  if(DEBUG){
    printf("DEBUG MODE: saving neighbor mesh\n");
    saveMesh(std::vector<Vector3f>(newPoints,newPoints+numPoints),neighbor_edges,"DEBUG_neighbors.obj");
    printf("DEBUG MODE: done.\n");
  }

  std::sort(neighbor_edges.begin(),neighbor_edges.end()); //sort min to max weight
  //substep: create MST rooted at prev idx
  printf("setting up mst...\n");
  std::vector<Edge> mst_edges;
  int* pointForests = (int*) malloc(numPoints*sizeof(int));//kruskal's alogrithm connects forests until done
#ifdef USE_OMP
  #pragma omp parallel for
#endif
  for(int i=0;i<numPoints;i++) pointForests[i] = i;
  for(unsigned int i=0;i<neighbor_edges.size();i++){ //going from min to max
    Edge curedge = neighbor_edges[i];
    int f1 = pointForests[curedge.v1];
    int f2 = pointForests[curedge.v2];
    if(f1!=f2){ //different forest
      mst_edges.push_back(curedge);
      int minNum = (f1<f2)? f1 : f2;
      for(int i=0;i<numPoints;i++){
        if(pointForests[i]==f1 || pointForests[i]==f2) pointForests[i] = minNum;
      }
    }
  }
  free(pointForests);
  
  if(DEBUG){
    printf("DEBUG MODE: saving MST mesh\n");
    saveMesh(std::vector<Vector3f>(newPoints,newPoints+numPoints),mst_edges,"DEBUG_MST.obj");
    printf("DEBUG MODE: done.\n");
  }

  //substep: create adjacency list for faster lookup
  printf("making adjacency list...\n");
  std::vector< std::vector<int> > adj(numPoints,std::vector<int>(0));
  for(unsigned int i=0;i<mst_edges.size();i++){
    Edge e = mst_edges[i];
    adj[e.v1].push_back(e.v2);
    adj[e.v2].push_back(e.v1);
  }
  
  //substep: use DFS from curIdx on mst to propogate normal directions
  printf("traversing mst to propagate normals...\n");
  bool* seen_mask = (bool*) calloc(numPoints,sizeof(bool));//calloc inits all to false
  std::vector<int> queue;
  queue.push_back(curidx);
  
  while(!queue.empty()){
    curidx = queue.back();
    queue.pop_back();
    seen_mask[curidx] = true; //processed
    for(unsigned int i=0;i<adj[curidx].size();i++){
      int nextidx = adj[curidx][i];
      if(!seen_mask[nextidx]){
        if(planes[curidx].normal.dot(planes[nextidx].normal)<0) planes[nextidx].normal*= -1.f; //flip dir
        queue.push_back(nextidx);
      }
    }
  }
  free(seen_mask);

  if(DEBUG){
    printf("DEBUG MODE: saving normals mesh\n");
    savePlanes(planes,numPoints,"DEBUG_normals.obj");
    printf("DEBUG MODE: done.\n");
  }

  //step 3: create a bounding box of the universe and split it into cubes
  printf("Approximating mesh now!\n");
  bbox system;
  for(int i=0;i<numPoints;i++) system.expand(points[i]);
  system.print();
  
  Vector3f universeSize = system.max-system.min;
  float sideLength = rho+delta; //cube side length
  //modify system size to fit int num of cubes in each dir
  printPoint(universeSize);
  float widthDif = int(std::ceil(universeSize(0)/sideLength))*sideLength-universeSize(0),
       heightDif = int(std::ceil(universeSize(1)/sideLength))*sideLength-universeSize(1),
        depthDif = int(std::ceil(universeSize(2)/sideLength))*sideLength-universeSize(2);
  Vector3f addon = Vector3f(widthDif,heightDif,depthDif);
  printPoint(addon);
  addon = addon/2.0; // add evenly to min and max
  system.max += addon;
  system.min -= addon;
  system.print();
  Vector3f numCubes = system.max-system.min;
  numCubes = numCubes/sideLength;
  if(numCubes(0)==0)numCubes(0)++;
  if(numCubes(1)==0)numCubes(1)++;
  if(numCubes(2)==0)numCubes(2)++;
  printf("Num cubes in each dir is:\n");
  printPoint(numCubes);

  //step 4: Approximate mesh based on differences between cubes
  bbox*** cubes = (bbox***) malloc(numCubes(0)*sizeof(bbox**));
  std::vector<Vector3f> newvertices;
  std::vector<Edge> edges;
  for(int i=0;i<numCubes(0);i++){
    cubes[i] = (bbox**) malloc(numCubes(1)*sizeof(bbox*));
    for(int j=0;j<numCubes(1);j++){
      cubes[i][j] = (bbox*) malloc(numCubes(2)*sizeof(bbox));
      for(int k=0;k<numCubes(2);k++){
        bbox cube = bbox(system.min+Vector3f(i*sideLength,j*sideLength,k*sideLength),
                         system.min+Vector3f((i+1)*sideLength,(j+1)*sideLength,(k+1)*sideLength));
        cubes[i][j][k] = cube;
        //printf("cube %d,%d,%d is:\n",i,j,k);cubes[i][j][k].print();
        Vector3f blb,blf,brb,brf,tlb,tlf,trb,trf; //[top/bottom][left/right][front/back] values at each corner
        float blbv,blfv,brbv,brfv,tlbv,tlfv,trbv,trfv; //actual vals at point
        blb=cube.min;
        blbv=getDist(blb,planes,numPoints);
        blf=cube.min+Vector3f(0,0,sideLength);
        blfv=getDist(blf,planes,numPoints);
        brb=cube.min+Vector3f(0,sideLength,0);
        brbv=getDist(brb,planes,numPoints);
        brf=cube.max-Vector3f(sideLength,0,0);
        brfv=getDist(brf,planes,numPoints);
        tlb=cube.min+Vector3f(sideLength,0,0);
        tlbv=getDist(tlb,planes,numPoints);
        tlf=cube.max-Vector3f(0,sideLength,0);
        tlfv=getDist(tlf,planes,numPoints);
        trb=cube.max-Vector3f(0,0,sideLength);
        trbv=getDist(trb,planes,numPoints);
        trf=cube.max;
        trfv=getDist(trf,planes,numPoints);
        for(int side=0;side<6;side++){//all six faces of cube
          Vector3f p[4]; float v[4];
          switch(side){ //grab points clockwise 
            case 0://top
              p[0]=tlb;p[1]=trb;p[2]=trf;p[3]=tlf;
              v[0]=tlbv;v[1]=trbv;v[2]=trfv;v[3]=tlfv;
              break;   
            case 1://bottom
              p[0]=blb;p[1]=brb;p[2]=brf;p[3]=blf;
              v[0]=blbv;v[1]=brbv;v[2]=brfv;v[3]=blfv;
              break; 
            case 2://left
              p[0]=tlb;p[1]=tlf;p[2]=blf;p[3]=blb;
              v[0]=tlbv;v[1]=tlfv;v[2]=blfv;v[3]=blbv;
              break;
            case 3://right
              p[0]=trf;p[1]=trb;p[2]=brb;p[3]=brf;
              v[0]=trfv;v[1]=trbv;v[2]=brbv;v[3]=brfv;
              break;
            case 4://front
              p[0]=tlf;p[1]=trf;p[2]=brf;p[3]=blf;
              v[0]=tlfv;v[1]=trfv;v[2]=brfv;v[3]=blfv;
              break;
            case 5://back
              p[0]=tlb;p[1]=trb;p[2]=brf;p[3]=blf;
              v[0]=tlbv;v[1]=trbv;v[2]=brfv;v[3]=blfv;
              break;
          }

          if(((v[0]<=0)==(v[1]<=0)) && ((v[1]<=0)==(v[2]<=0)) && ((v[2]<=0)==(v[3]<=0))) continue; //4 match means not part of mesh
          else if(((v[0]<=0)==(v[1]<=0) && (v[1]<=0)==(v[2]<=0)) ||
                  ((v[0]<=0)==(v[1]<=0) && (v[1]<=0)==(v[3]<=0)) ||
                  ((v[0]<=0)==(v[2]<=0) && (v[2]<=0)==(v[3]<=0)) ||
                  ((v[1]<=0)==(v[2]<=0) && (v[2]<=0)==(v[3]<=0))){//three match
            for(int idx=0;idx<4;idx++){
              float vm = v[idx];Vector3f pm = p[idx];
              float vp = v[(idx-1)%4];Vector3f pp = p[(idx-1)%4];
              float vn = v[(idx+1)%4];Vector3f pn = p[(idx+1)%4];
              if((vm<=0)==(vp<=0) || (vm<=0)==(vn<=0)) continue; //only consider odd man out
              float frac1 = fabs(vm)/(fabs(vm)+fabs(vp)),
                    frac2 = fabs(vm)/(fabs(vm)+fabs(vn));
              Vector3f newP1 = (pm*(1-frac1))+(pp*frac1);
              Vector3f newP2 = (pm*(1-frac2))+(pn*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(Edge(newvertices.size()-1,newvertices.size()-2));
            }
          }else{//two match
            if((v[0]<=0)==(v[1]<=0)){ //horiz line
              float frac1 = fabs(v[0])/(fabs(v[0])+fabs(v[3])),
                    frac2 = fabs(v[1])/(fabs(v[1])+fabs(v[2]));
              Vector3f newP1 = (p[0]*(1-frac1))+(p[3]*frac1);
              Vector3f newP2 = (p[1]*(1-frac2))+(p[2]*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(Edge(newvertices.size()-1,newvertices.size()-2));
            }else if((v[0]<=0)==(v[3]<=0)){ //vert line
              float frac1 = fabs(v[0])/(fabs(v[0])+fabs(v[1])),
                    frac2 = fabs(v[3])/(fabs(v[3])+fabs(v[2]));
              Vector3f newP1 = (p[0]*(1-frac1))+(p[1]*frac1);
              Vector3f newP2 = (p[3]*(1-frac2))+(p[2]*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(Edge(newvertices.size()-1,newvertices.size()-2));
            }else{ //double diagonal
              float frac1 = fabs(v[0])/(fabs(v[0])+fabs(v[1])),
                    frac2 = fabs(v[0])/(fabs(v[0])+fabs(v[3])),
                    frac3 = fabs(v[2])/(fabs(v[2])+fabs(v[1])),
                    frac4 = fabs(v[2])/(fabs(v[2])+fabs(v[3]));
              Vector3f newP1 = (p[0]*(1-frac1))+(p[1]*frac1);
              Vector3f newP2 = (p[0]*(1-frac2))+(p[3]*frac2);
              Vector3f newP3 = (p[2]*(1-frac3))+(p[1]*frac3);
              Vector3f newP4 = (p[2]*(1-frac4))+(p[3]*frac4);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              newvertices.push_back(newP3);
              newvertices.push_back(newP4);

              //edge dir depends on center val;
              Vector3f pcent = 0.25*(newP1+newP2+newP3+newP4);
              float vcent = getDist(pcent,planes,numPoints);
              if((vcent<=0)==(v[0]<=0)){
                edges.push_back(Edge(newvertices.size()-4,newvertices.size()-2));
                edges.push_back(Edge(newvertices.size()-1,newvertices.size()-3));
              }else{
                edges.push_back(Edge(newvertices.size()-3,newvertices.size()-4));
                edges.push_back(Edge(newvertices.size()-1,newvertices.size()-2));
              }
            }
          }
        }
      }
    }
  }
  free(planes);

  //get rid of duplicate vertices
  int* mapping = (int*) malloc(sizeof(int)*newvertices.size());
  for(unsigned int i=0;i<newvertices.size();i++){
    Vector3f v1 = newvertices[i];
    bool unique = true;
    for(unsigned int j=0;j<finalVertices.size();j++){
      if(v1==finalVertices[j]){
        unique = false;
        mapping[i] = j;
        break;
      }
    }
    if(unique){
      mapping[i] = finalVertices.size();
      finalVertices.push_back(v1);
    }
  }
  
  //map edges onto unique vertices, delete duplicate edges
  bool** seenMat = (bool**) malloc(sizeof(bool*)*finalVertices.size());
  for(unsigned int i=0;i<finalVertices.size(); i++) seenMat[i] = (bool*) calloc(finalVertices.size(),sizeof(bool));

  for(unsigned int i=0;i<edges.size();i++){
    int v1=mapping[edges[i].v1],
        v2=mapping[edges[i].v2];
    if(!seenMat[v1][v2]){
      seenMat[v1][v2]=true;
      finalEdges.push_back(Edge(v1,v2));
    }
  }
  if(DEBUG){
    printf("DEBUG MODE: saving approximated mesh\n");
    saveMesh(finalVertices,finalEdges,"DEBUG_approximate.obj");
    printf("DEBUG MODE: done.\n");
  }
}
