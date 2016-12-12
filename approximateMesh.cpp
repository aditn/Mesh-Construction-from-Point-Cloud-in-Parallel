/* needed for malloc, printf */
#include <stdio.h>
#include <stdlib.h>

/* needed for vectors */
#include <vector>
#include <algorithm> //for sorting vectors

/* fabs, square, sqrt */
#include <cmath>

#ifdef USE_OMP
/* needed for openMP! */
#include <omp.h>
#endif

/* data structures */
#include "structures.h"
#include "tangentPlane.h"
#include "parseOBJ.h"
#include "approximateMesh.h"
#include "ourTimer.h"
#include "constants.h"

/* linear algebra library */
#include <Eigen/Dense>

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

inline bool sameSide(float val1,float val2){
  //assumes val is the dist
  //returns true if have the same sign
  return (val1>=0) == (val2>=0);
}

inline Vector3f interpolate(Vector3f pA,Vector3f pB,float vA,float vB){
  float frac = fabs(vA)/(fabs(vA)+fabs(vB));
  return pA+frac*(pB-pA);
}

void approximateMesh(Vector3f* points, int numPoints,std::vector<Vector3f>& finalVertices,std::vector<Edge>& finalEdges){
  //step 1: create a plane for each point based off of its neighbors
  printf("getting planes...\n");
#ifdef USE_OMP
  printf("using OMP with %d threads\n",numThreads);
#endif
  Plane* planes = computeTangentPlanes(points,numPoints);
  Vector3f* newPoints = (Vector3f*) malloc(sizeof(Vector3f)*numPoints); // Plane centroids
#ifdef USE_OMP
  #pragma omp parallel for
#endif
  for(int i=0;i<numPoints;i++) newPoints[i] = planes[i].center;
  printf("time %.4fs\n",timeSince());
  
  //step 2: propogate normal directions of every plane
  //substep: find plane centroid with largest z coordinate
  printf("getting max z...\n");
  int curmax = 0;
  int curidx = -1;
  /* Starting point is point with largest z coordinate */
  for(int i=0;i<numPoints;i++){
    if(curidx==-1 || newPoints[i](2)>curmax){
      curmax = newPoints[i](2);
      curidx = i;
    }
  }
  planes[curidx].normal = Vector3f(0,0,(curmax>=0)? -1: 1); //TODO: why?

  //substep: create graph of neighbors with edge weights (MST)
  printf("getting neighbor graph...\n");
#ifdef USE_OMP
  std::vector<int>* neighbors = (std::vector<int>*) malloc(sizeof(std::vector<Edge>)*numPoints);
  #pragma omp parallel for
  for(int i=0;i<numPoints;i++){
    neighbors[i] = getNearest(newPoints,numPoints,i,K);
  }
  std::vector<Edge> neighbor_edges;
  for(int i=0;i<numPoints;i++){
    for(unsigned int j=0;j<neighbors[i].size();j++){
      int idx = neighbors[i][j];
      neighbor_edges.push_back(Edge(i,idx,1-fabs(planes[i].normal.dot(planes[idx].normal))));
    }
  }
  free(neighbors);
#else
  std::vector<Edge> neighbor_edges;
  for(int i=0;i<numPoints;i++){
    std::vector<int> neighbors = getNearest(newPoints,numPoints,i,K);
    for(unsigned int j=0;j<neighbors.size();j++){
      int idx = neighbors[j];
      neighbor_edges.push_back(Edge(i,idx,1-fabs(planes[i].normal.dot(planes[idx].normal))));
    }
  }
#endif
  printf("time %.4fs\n",timeSince());
  
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
  printf("time %.4fs\n",timeSince());
  
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
  printf("time %.4fs\n",timeSince());

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
  float sideLength = ro+delta; //cube side length
  //modify system size to fit int num of cubes in each dir
  float widthDif = int(std::ceil(universeSize(0)/sideLength))*sideLength-universeSize(0),
       heightDif = int(std::ceil(universeSize(1)/sideLength))*sideLength-universeSize(1),
        depthDif = int(std::ceil(universeSize(2)/sideLength))*sideLength-universeSize(2);
  Vector3f addon = Vector3f(widthDif,heightDif,depthDif);
  addon = addon*0.5; // add evenly to min and max
  system.max += addon;
  system.min -= addon;
  Vector3f numCubes = system.max-system.min;
  numCubes = numCubes/sideLength;
  if(!numCubes(0)) numCubes(0)++;
  if(!numCubes(1)) numCubes(1)++;
  if(!numCubes(2)) numCubes(2)++;
  printf("Num cubes in each dir is:\n");
  printPoint(numCubes);

  //step 4: Approximate mesh based on differences between cubes

#ifdef USE_OMP
  std::vector<Vector3f> newvertices;
  std::vector<Edge> edges;

  const int cubecount = numCubes(0)*numCubes(1)*numCubes(2);
  const int xyCount = numCubes(0)*numCubes(1);
  const int xcount = numCubes(0);

  const int MAX_COUNT = 24; //max # verts/edges per cube
  Vector3f* percube_vertices = (Vector3f*) malloc(sizeof(Vector3f)*MAX_COUNT*cubecount);
  int* percube_vertex_counts = (int*) malloc(sizeof(int)*cubecount);
  Edge* percube_edges = (Edge*) malloc(sizeof(Vector3f)*MAX_COUNT*cubecount);
  int* percube_edge_counts = (int*) malloc(sizeof(int)*cubecount);

  #pragma omp parallel for
  for(int cubeIdx=0;cubeIdx<cubecount;cubeIdx++){
    int intermediate = cubeIdx % xyCount;
    int k = cubeIdx / xyCount;
    int j = intermediate / xcount;
    int i = intermediate % xcount;
    bbox cube = bbox(system.min+Vector3f(i*sideLength,j*sideLength,k*sideLength),
                     system.min+Vector3f((i+1)*sideLength,(j+1)*sideLength,(k+1)*sideLength));
    Vector3f blb,blf,brb,brf,tlb,tlf,trb,trf; //[top/bottom][left/right][front/back] values at each corner
    float blbv,blfv,brbv,brfv,tlbv,tlfv,trbv,trfv; //actual vals at point
    blb=cube.min;
    blbv=getDist(blb,planes,numPoints);
    blf=cube.min+Vector3f(0,0,sideLength);
    blfv=getDist(blf,planes,numPoints);
    brb=cube.min+Vector3f(sideLength,0,0);
    brbv=getDist(brb,planes,numPoints);
    brf=cube.max-Vector3f(0,sideLength,0);
    brfv=getDist(brf,planes,numPoints);
    tlb=cube.min+Vector3f(0,sideLength,0);
    tlbv=getDist(tlb,planes,numPoints);
    tlf=cube.max-Vector3f(sideLength,0,0);
    tlfv=getDist(tlf,planes,numPoints);
    trb=cube.max-Vector3f(0,0,sideLength);
    trbv=getDist(trb,planes,numPoints);
    trf=cube.max;
    trfv=getDist(trf,planes,numPoints);

    int vertex_count = 0;
    int edge_count = 0;
    const int base = MAX_COUNT*cubeIdx;

    for(int side=0;side<6;side++){//all six faces of cube
      Vector3f p[4]; float v[4];
      switch(side){ //grab points clockwise 
       case 0://top
         p[0]=tlb; p[1]=trb; p[2]=trf; p[3]=tlf;
          v[0]=tlbv;v[1]=trbv;v[2]=trfv;v[3]=tlfv;
          break;   
        case 1://bottom
          p[0]=blb; p[1]=brb; p[2]=brf; p[3]=blf;
          v[0]=blbv;v[1]=brbv;v[2]=brfv;v[3]=blfv;
          break; 
        case 2://left
          p[0]=tlb; p[1]=tlf; p[2]=blf; p[3]=blb;
          v[0]=tlbv;v[1]=tlfv;v[2]=blfv;v[3]=blbv;
          break;
        case 3://right
          p[0]=trb; p[1]=trf; p[2]=brf; p[3]=brb;
          v[0]=trbv;v[1]=trfv;v[2]=brfv;v[3]=brbv;
          break;
        case 4://front
          p[0]=tlf; p[1]=trf; p[2]=brf; p[3]=blf;
          v[0]=tlfv;v[1]=trfv;v[2]=brfv;v[3]=blfv;
          break;
        case 5://back
          p[0]=tlb; p[1]=trb; p[2]=brb; p[3]=blb;
          v[0]=tlbv;v[1]=trbv;v[2]=brbv;v[3]=blbv;
          break;
      }

      if(sameSide(v[0],v[1]) && sameSide(v[1],v[2]) && sameSide(v[2],v[3])) continue; //no edges needed
      else if( (sameSide(v[0],v[1]) && sameSide(v[1],v[2])) ||
               (sameSide(v[0],v[1]) && sameSide(v[1],v[3])) ||
               (sameSide(v[0],v[2]) && sameSide(v[2],v[3])) ||
               (sameSide(v[1],v[2]) && sameSide(v[2],v[3]))){//three match
        for(int idx=0;idx<4;idx++){
          float vm = v[idx];Vector3f pm = p[idx];
          float vp = v[(idx+3)%4];Vector3f pp = p[(idx+3)%4];
          float vn = v[(idx+1)%4];Vector3f pn = p[(idx+1)%4];
          if(sameSide(vm,vp) || sameSide(vm,vn)) continue; //only consider odd man out
          percube_vertices[base+(vertex_count++)] = interpolate(pm,pp,vm,vp);
          percube_vertices[base+(vertex_count++)] = interpolate(pm,pn,vm,vn);
          percube_edges[base+(edge_count++)] = Edge(vertex_count-1,vertex_count-2);
        }
      }else{//two match
        if(sameSide(v[0],v[1])){ //horiz line
          percube_vertices[base+(vertex_count++)]=interpolate(p[0],p[3],v[0],v[3]);
          percube_vertices[base+(vertex_count++)]=interpolate(p[1],p[2],v[1],v[2]);
          percube_edges[base+(edge_count++)]=Edge(vertex_count-1,vertex_count-2);
        }else if(sameSide(v[0],v[3])){ //vert line
          percube_vertices[base+(vertex_count++)]=interpolate(p[0],p[1],v[0],v[1]);
          percube_vertices[base+(vertex_count++)]=interpolate(p[3],p[2],v[3],v[2]);
          percube_edges[base+(edge_count++)]=Edge(vertex_count-1,vertex_count-2);
        }else{ //double diagonal
          Vector3f newP1 = interpolate(p[0],p[1],v[0],v[1]);
          Vector3f newP2 = interpolate(p[0],p[3],v[0],v[3]);
          Vector3f newP3 = interpolate(p[2],p[1],v[2],v[1]);
          Vector3f newP4 = interpolate(p[2],p[3],v[2],v[3]);
          percube_vertices[base+(vertex_count++)]=newP1;
          percube_vertices[base+(vertex_count++)]=newP2;
          percube_vertices[base+(vertex_count++)]=newP3;
          percube_vertices[base+(vertex_count++)]=newP4;
  
          //edge dir depends on center val;
          Vector3f pcenter = 0.25*(newP1+newP2+newP3+newP4);
          if(sameSide(getDist(pcenter,planes,numPoints),v[0])){
            percube_edges[base+(edge_count++)]=Edge(vertex_count-4,vertex_count-2);
            percube_edges[base+(edge_count++)]=Edge(vertex_count-1,vertex_count-3);
          }else{
            percube_edges[base+(edge_count++)]=Edge(vertex_count-3,vertex_count-4);
            percube_edges[base+(edge_count++)]=Edge(vertex_count-1,vertex_count-2);
          }
        }
      }
    }
    percube_vertex_counts[cubeIdx]=vertex_count;
    percube_edge_counts[cubeIdx]=edge_count;
  }

  //coalesce edge/vertex 'packets'
  int vertexOffset = 0;
  for(int i=0;i<cubecount;i++){
    const int base = MAX_COUNT*i;
    for(int j=0;j<percube_vertex_counts[i];j++){
      newvertices.push_back(percube_vertices[base+j]);
    }
    for(int j=0;j<percube_edge_counts[i];j++){
      Edge cur_edge = percube_edges[base+j];
      cur_edge.v1 += vertexOffset;
      cur_edge.v2 += vertexOffset;
      edges.push_back(cur_edge);
    }
    vertexOffset += percube_vertex_counts[i];
  }

#else
  std::vector<Vector3f> newvertices;
  std::vector<Edge> edges;
  for(int i=0;i<numCubes(0);i++){
    for(int j=0;j<numCubes(1);j++){
      for(int k=0;k<numCubes(2);k++){
        bbox cube = bbox(system.min+Vector3f(i*sideLength,j*sideLength,k*sideLength),
                         system.min+Vector3f((i+1)*sideLength,(j+1)*sideLength,(k+1)*sideLength));
        //printf("cube %d,%d,%d is:\n",i,j,k);cubes[i][j][k].print();
        Vector3f blb,blf,brb,brf,tlb,tlf,trb,trf; //[top/bottom][left/right][front/back] values at each corner
        float blbv,blfv,brbv,brfv,tlbv,tlfv,trbv,trfv; //actual vals at point
        blb=cube.min;
        blbv=getDist(blb,planes,numPoints);
        blf=cube.min+Vector3f(0,0,sideLength);
        blfv=getDist(blf,planes,numPoints);
        brb=cube.min+Vector3f(sideLength,0,0);
        brbv=getDist(brb,planes,numPoints);
        brf=cube.max-Vector3f(0,sideLength,0);
        brfv=getDist(brf,planes,numPoints);
        tlb=cube.min+Vector3f(0,sideLength,0);
        tlbv=getDist(tlb,planes,numPoints);
        tlf=cube.max-Vector3f(sideLength,0,0);
        tlfv=getDist(tlf,planes,numPoints);

        trb=cube.max-Vector3f(0,0,sideLength);
        trbv=getDist(trb,planes,numPoints);
        trf=cube.max;
        trfv=getDist(trf,planes,numPoints);
        for(int side=0;side<6;side++){//all six faces of cube
          Vector3f p[4]; float v[4];
          switch(side){ //grab points clockwise 
            case 0://top
              p[0]=tlb; p[1]=trb; p[2]=trf; p[3]=tlf;
              v[0]=tlbv;v[1]=trbv;v[2]=trfv;v[3]=tlfv;
              break;   
            case 1://bottom
              p[0]=blb; p[1]=brb; p[2]=brf; p[3]=blf;
              v[0]=blbv;v[1]=brbv;v[2]=brfv;v[3]=blfv;
              break; 
            case 2://left
              p[0]=tlb; p[1]=tlf; p[2]=blf; p[3]=blb;
              v[0]=tlbv;v[1]=tlfv;v[2]=blfv;v[3]=blbv;
              break;
            case 3://right
              p[0]=trb; p[1]=trf; p[2]=brf; p[3]=brb;
              v[0]=trbv;v[1]=trfv;v[2]=brfv;v[3]=brbv;
              break;
            case 4://front
              p[0]=tlf; p[1]=trf; p[2]=brf; p[3]=blf;
              v[0]=tlfv;v[1]=trfv;v[2]=brfv;v[3]=blfv;
              break;
            case 5://back
              p[0]=tlb; p[1]=trb; p[2]=brb; p[3]=blb;
              v[0]=tlbv;v[1]=trbv;v[2]=brbv;v[3]=blbv;
              break;
          }

          if(sameSide(v[0],v[1]) && sameSide(v[1],v[2]) && sameSide(v[2],v[3])) continue; //no edges needed
          else if( (sameSide(v[0],v[1]) && sameSide(v[1],v[2])) ||
                   (sameSide(v[0],v[1]) && sameSide(v[1],v[3])) ||
                   (sameSide(v[0],v[2]) && sameSide(v[2],v[3])) ||
                   (sameSide(v[1],v[2]) && sameSide(v[2],v[3]))){//three match
            for(int idx=0;idx<4;idx++){
              float vm = v[idx];Vector3f pm = p[idx];
              float vp = v[(idx+3)%4];Vector3f pp = p[(idx+3)%4];
              float vn = v[(idx+1)%4];Vector3f pn = p[(idx+1)%4];
              if(sameSide(vm,vp) || sameSide(vm,vn)) continue; //only consider odd man out
              newvertices.push_back(interpolate(pm,pp,vm,vp));
              newvertices.push_back(interpolate(pm,pn,vm,vn));
              edges.push_back(Edge(newvertices.size()-1,newvertices.size()-2));
            }
          }else{//two match
            if(sameSide(v[0],v[1])){ //horiz line
              newvertices.push_back(interpolate(p[0],p[3],v[0],v[3]));
              newvertices.push_back(interpolate(p[1],p[2],v[1],v[2]));
              edges.push_back(Edge(newvertices.size()-1,newvertices.size()-2));
            }else if(sameSide(v[0],v[3])){ //vert line
              newvertices.push_back(interpolate(p[0],p[1],v[0],v[1]));
              newvertices.push_back(interpolate(p[3],p[2],v[3],v[2]));
              edges.push_back(Edge(newvertices.size()-1,newvertices.size()-2));
            }else{ //double diagonal
              Vector3f newP1 = interpolate(p[0],p[1],v[0],v[1]);
              Vector3f newP2 = interpolate(p[0],p[3],v[0],v[3]);
              Vector3f newP3 = interpolate(p[2],p[1],v[2],v[1]);
              Vector3f newP4 = interpolate(p[2],p[3],v[2],v[3]);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              newvertices.push_back(newP3);
              newvertices.push_back(newP4);

              //edge dir depends on center val;
              Vector3f pcenter = 0.25*(newP1+newP2+newP3+newP4);
              if(sameSide(getDist(pcenter,planes,numPoints),v[0])){
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
#endif
  free(planes);
  printf("time %.4fs\n",timeSince());

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
