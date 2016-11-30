/* needed for malloc, printf */
#include <stdio.h>
#include <stdlib.h>

/* needed for vectors */
#include <vector>

/* fabs, square, sqrt */
#include <cmath>

/* data structures */
#include "structures.h"
#include "tangentPlane.h"
#include "approximateMesh.h"

void approximateMesh(V3* points, int numPoints, float rho, float delta,std::vector<V3>& finalVertices,std::vector<Edge>& finalEdges){
  //step 1: create a plane for each point based off of its neighbors
  Plane* planes = computeTangentPlanes(points,numPoints,rho,delta);

  //step 2: propogate normal directions of every plane
  
  //step 3: create a bounding box of the universe and split it into cubes
  bbox system;
  for(int i=0;i<numPoints;i++) system.expand(points[i]);
  system.print();
  
  V3 universeSize = system.max-system.min;
  float sideLength = rho+delta; //cube side length
  //modify system size to fit int num of cubes in each dir
  printPoint(universeSize);
  float widthDif = universeSize.x-int(std::ceil(universeSize.x/sideLength))*sideLength,
     heightDif = universeSize.y-int(std::ceil(universeSize.y/sideLength))*sideLength,
      depthDif = universeSize.z-int(std::ceil(universeSize.z/sideLength))*sideLength;
  V3 addon = V3(-widthDif,-heightDif,-depthDif);
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

  //step 2e: Approximate mesh based on differences between cubes
  bbox*** cubes = (bbox***) malloc(numCubes.x*sizeof(bbox**));
  std::vector<V3> newvertices;
  std::vector<E> edges;
  for(int i=0;i<numCubes.x;i++){
    cubes[i] = (bbox**) malloc(numCubes.y*sizeof(bbox*));
    for(int j=0;j<numCubes.y;j++){
      cubes[i][j] = (bbox*) malloc(numCubes.z*sizeof(bbox));
      for(int k=0;k<numCubes.z;k++){
        bbox cube = bbox(system.min+V3(i*sideLength,j*sideLength,k*sideLength),sideLength,sideLength,sideLength);
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
        for(int side=0;side<6;side++){//every face of cube
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

          if((v[0]<=0)==(v[1]<=0)==(v[2]<=0)==(v[3]<=0)) continue; //4 match
          else if(((v[0]<=0)==(v[1]<=0) && (v[1]<=0)==(v[2]<=0)) ||
                  ((v[0]<=0)==(v[1]<=0) && (v[1]<=0)==(v[3]<=0)) ||
                  ((v[0]<=0)==(v[2]<=0) && (v[2]<=0)==(v[3]<=0)) ||
                  ((v[1]<=0)==(v[2]<=0) && (v[2]<=0)==(v[3]<=0))){//three match
            for(int idx=0;idx<4;idx++){
              float vm = v[idx];V3 pm = p[idx];
              float vp = v[(idx-1)%4];V3 pp = p[(idx-1)%4];
              float vn = v[(idx+1)%4];V3 pn = p[(idx+1)%4];
              if((vm<=0)==(vp<=0) || (vm<=0)==(vn<=0)) continue; //only consider odd man out
              float frac1 = fabs(vm)/fabs(vm+vp),
                    frac2 = fabs(vm)/fabs(vm+vn);
              V3 newP1 = pm+(pp*frac1);
              V3 newP2 = pm+(pn*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(E(newP1,newP2));
            }
          }else{//two match
            if((v[0]<=0)==(v[1]<=0)){ //horiz line
              float frac1 = fabs(v[0])/fabs(v[0]+v[3]),
                    frac2 = fabs(v[1])/fabs(v[1]+v[2]);
              V3 newP1 = p[0]+(p[3]*frac1);
              V3 newP2 = p[1]+(p[2]*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(E(newP1,newP2));
            }else if((v[0]<=0)==(v[3]<=0)){ //vert line
              float frac1 = fabs(v[0])/fabs(v[0]+v[1]),
                    frac2 = fabs(v[3])/fabs(v[3]+v[2]);
              V3 newP1 = p[0]+(p[1]*frac1);
              V3 newP2 = p[3]+(p[2]*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(E(newP1,newP2));
            }else{ //double diagonal
              float frac1 = fabs(v[0])/fabs(v[0]+v[1]),
                    frac2 = fabs(v[0])/fabs(v[0]+v[3]),
                    frac3 = fabs(v[2])/fabs(v[2]+v[1]),
                    frac4 = fabs(v[2])/fabs(v[2]+v[3]);
              V3 newP1 = p[0]+(p[1]*frac1);
              V3 newP2 = p[0]+(p[3]*frac2);
              V3 newP3 = p[2]+(p[1]*frac3);
              V3 newP4 = p[2]+(p[3]*frac4);
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
