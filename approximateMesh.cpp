/* needed for malloc, printf */
#include <stdio.h>
#include <stdlib.h>

/* needed for vectors */
#include <vector>
#include <algorithm> //for sorting vectors

/* fabs, square, sqrt */
#include <cmath>

/* data structures */
#include "structures.h"
#include "tangentPlane.h"
#include "approximateMesh.h"

#define K 5 //num neighbors for MST Propogation

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

std::vector<int> getNearest(V3* points, int numPoints, int idx, int numNeighbors){
  std::vector<int> neighbors;
  if(numPoints<=numNeighbors){
    for(int i=0;i<numPoints;i++) neighbors.push_back(i);
    return neighbors;
  }else{
    int* nearestIdxs = (int*) malloc(sizeof(int)*numNeighbors);
    float* nearestVals = (float*) malloc(sizeof(float)*numNeighbors);
    for(int i=0;i<numNeighbors;i++) nearestVals[i] = INFINITY;

    V3 curPoint = points[idx];
    for(int i=0;i<numPoints;i++){
      if(i==idx) continue;
      insMin(nearestIdxs,nearestVals,numNeighbors,i,curPoint.dist(points[i]));
    }
    for(int i=0;i<numNeighbors;i++) neighbors.push_back(nearestIdxs[i]);
    free(nearestIdxs);
    free(nearestVals);
    return neighbors;
  }
}

void approximateMesh(V3* points, int numPoints, float rho, float delta,std::vector<V3>& finalVertices,std::vector<Edge>& finalEdges){
  //step 1: create a plane for each point based off of its neighbors
  printf("getting planes...\n");
  Plane* planes = computeTangentPlanes(points,numPoints,rho,delta);
  V3* newPoints = (V3*) malloc(sizeof(V3)*numPoints);
  for(int i=0;i<numPoints;i++) newPoints[i] = planes[i].center;
  //step 2: propogate normal directions of every plane
  //substep: find plane centroid with largest z coordinate
  printf("getting max z...\n");
  int curmax = 0;
  int curidx = -1;
  for(int i=0;i<numPoints;i++){
    if(curidx==-1 || newPoints[i].z>curmax){
      curmax = newPoints[i].z;
      curidx = i;
    }
  }
  planes[curidx].normal = V3(0,0,(curmax>=0)? -1: 1); //TODO: why?

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

  std::sort(neighbor_edges.begin(),neighbor_edges.end()); //sort min to max weight
  //substep: create MST rooted at prev idx
  printf("setting up mst...\n");
  std::vector<Edge> mst_edges;
  int* pointForests = (int*) malloc(numPoints*sizeof(int));//kruskal's alogrithm connects forests until done
  for(int i=0;i<numPoints;i++) pointForests[i] = i;
  for(unsigned int i=0;i<neighbor_edges.size();i++){ //going from min to max
    Edge curedge = neighbor_edges[i];
    int f1 = pointForests[curedge.v1];
    int f2 = pointForests[curedge.v2];
    if(f1!=f2){ //different forest
      mst_edges.push_back(curedge);
      if(f1<f2) for(int i=0;i<numPoints;i++) if(pointForests[i]==f2) pointForests[i]=f1;
      else      for(int i=0;i<numPoints;i++) if(pointForests[i]==f1) pointForests[i]=f2;
    }
  }
  free(pointForests);

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

  //step 3: create a bounding box of the universe and split it into cubes
  printf("Approximating mesh now!\n");
  bbox system;
  for(int i=0;i<numPoints;i++) system.expand(points[i]);
  system.print();
  
  V3 universeSize = system.max-system.min;
  float sideLength = rho+delta; //cube side length
  //modify system size to fit int num of cubes in each dir
  printPoint(universeSize);
  float widthDif = int(std::ceil(universeSize.x/sideLength))*sideLength-universeSize.x,
       heightDif = int(std::ceil(universeSize.y/sideLength))*sideLength-universeSize.y,
        depthDif = int(std::ceil(universeSize.z/sideLength))*sideLength-universeSize.z;
  V3 addon = V3(widthDif,heightDif,depthDif);
  printPoint(addon);
  addon.scale(0.5); // add evenly to min and max
  system.max += addon;
  system.min -= addon;
  system.print();
  V3 numCubes = system.max-system.min;
  numCubes.scale(1.0/sideLength);
  if(numCubes.x==0)numCubes.x++;
  if(numCubes.y==0)numCubes.y++;
  if(numCubes.z==0)numCubes.z++;
  printf("Num cubes in each dir is:\n");
  printPoint(numCubes);

  //step 4: Approximate mesh based on differences between cubes
  bbox*** cubes = (bbox***) malloc(numCubes.x*sizeof(bbox**));
  std::vector<V3> newvertices;
  std::vector<E> edges;
  for(int i=0;i<numCubes.x;i++){
    cubes[i] = (bbox**) malloc(numCubes.y*sizeof(bbox*));
    for(int j=0;j<numCubes.y;j++){
      cubes[i][j] = (bbox*) malloc(numCubes.z*sizeof(bbox));
      for(int k=0;k<numCubes.z;k++){
        bbox cube = bbox(system.min+V3(i*sideLength,j*sideLength,k*sideLength),
                         system.min+V3((i+1)*sideLength,(j+1)*sideLength,(k+1)*sideLength));
        cubes[i][j][k] = cube;
        //printf("cube %d,%d,%d is:\n",i,j,k);cubes[i][j][k].print();
        V3 blb,blf,brb,brf,tlb,tlf,trb,trf; //[top/bottom][left/right][front/back] values at each corner
        float blbv,blfv,brbv,brfv,tlbv,tlfv,trbv,trfv; //actual vals at point
        blb=cube.min;
        blbv=getDist(blb,planes,numPoints);
        blf=cube.min+V3(0,0,sideLength);
        blfv=getDist(blf,planes,numPoints);
        brb=cube.min+V3(0,sideLength,0);
        brbv=getDist(brb,planes,numPoints);
        brf=cube.max-V3(sideLength,0,0);
        brfv=getDist(brf,planes,numPoints);
        tlb=cube.min+V3(sideLength,0,0);
        tlbv=getDist(tlb,planes,numPoints);
        tlf=cube.max-V3(0,sideLength,0);
        tlfv=getDist(tlf,planes,numPoints);
        trb=cube.max-V3(0,0,sideLength);
        trbv=getDist(trb,planes,numPoints);
        trf=cube.max;
        trfv=getDist(trf,planes,numPoints);
        for(int side=0;side<6;side++){//all six faces of cube
          V3 p[4]; float v[4];
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
              float vm = v[idx];V3 pm = p[idx];
              float vp = v[(idx-1)%4];V3 pp = p[(idx-1)%4];
              float vn = v[(idx+1)%4];V3 pn = p[(idx+1)%4];
              if((vm<=0)==(vp<=0) || (vm<=0)==(vn<=0)) continue; //only consider odd man out
              float frac1 = fabs(vm)/(fabs(vm)+fabs(vp)),
                    frac2 = fabs(vm)/(fabs(vm)+fabs(vn));
              V3 newP1 = (pm*(1-frac1))+(pp*frac1);
              V3 newP2 = (pm*(1-frac2))+(pn*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(E(newP1,newP2));
            }
          }else{//two match
            if((v[0]<=0)==(v[1]<=0)){ //horiz line
              float frac1 = fabs(v[0])/(fabs(v[0])+fabs(v[3])),
                    frac2 = fabs(v[1])/(fabs(v[1])+fabs(v[2]));
              V3 newP1 = (p[0]*(1-frac1))+(p[3]*frac1);
              V3 newP2 = (p[1]*(1-frac2))+(p[2]*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(E(newP1,newP2));
            }else if((v[0]<=0)==(v[3]<=0)){ //vert line
              float frac1 = fabs(v[0])/(fabs(v[0])+fabs(v[1])),
                    frac2 = fabs(v[3])/(fabs(v[3])+fabs(v[2]));
              V3 newP1 = (p[0]*(1-frac1))+(p[1]*frac1);
              V3 newP2 = (p[3]*(1-frac2))+(p[2]*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(E(newP1,newP2));
            }else{ //double diagonal
              float frac1 = fabs(v[0])/(fabs(v[0])+fabs(v[1])),
                    frac2 = fabs(v[0])/(fabs(v[0])+fabs(v[3])),
                    frac3 = fabs(v[2])/(fabs(v[2])+fabs(v[1])),
                    frac4 = fabs(v[2])/(fabs(v[2])+fabs(v[3]));
              V3 newP1 = (p[0]*(1-frac1))+(p[1]*frac1);
              V3 newP2 = (p[0]*(1-frac2))+(p[3]*frac2);
              V3 newP3 = (p[2]*(1-frac3))+(p[1]*frac3);
              V3 newP4 = (p[2]*(1-frac4))+(p[3]*frac4);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              newvertices.push_back(newP3);
              newvertices.push_back(newP4);

              //edge dir depends on center val;
              V3 pcent = 0.25*(newP1+newP2+newP3+newP4);
              float vcent = getDist(pcent,planes,numPoints);
              if((vcent<=0)==(v[0]<=0)){
                edges.push_back(E(newP1,newP3));
                edges.push_back(E(newP2,newP4));
              }else{
                edges.push_back(E(newP1,newP2));
                edges.push_back(E(newP3,newP4));
              }
            }
          }
        }
      }
    }
  }
  free(planes);

  //get rid of duplicate vertices
  for(unsigned int i=0;i<newvertices.size();i++){
    V3 v1 = newvertices[i];
    bool unique = true;
    for(unsigned int j=0;j<finalVertices.size();j++){
      if(v1==finalVertices[j]){
        unique = false;
        break;
      }
    }
    if(unique) finalVertices.push_back(v1);
  }
  
  //go from edges as two vertices to edges as two indices into vertex buffer
  for(unsigned int i=0;i<edges.size();i++){
    V3 v1=edges[i].v1,
       v2=edges[i].v2;
    int v1x=-1,
        v2x=-1;
    for(unsigned int j=0;j<finalVertices.size();j++){
      if(v1x<0 && v1==finalVertices[j]) v1x=j;
      if(v2x<0 && v2==finalVertices[j]) v2x=j;
      if(v1x>=0 && v2x>=0) break;
    }
    if(v1x==-1 || v2x==-1) printf("Error. Something went wrong\n"); //couldn't find vertex in v buffer
    finalEdges.push_back(Edge(v1x,v2x));
  }
}
