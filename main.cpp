#define INPUT_FILE "inputs/sphere1000.obj"

/* needed for file processing */
#include <iostream>
#include <fstream>

/* needed for vectors */
#include <vector>
#include <algorithm> //for vector sorting

/* needed for malloc and printf and scanf and what not */
#include <stdio.h>
#include <stdlib.h>

/* sqrt and square */
#include <cmath>

/* data structures */
#include "structures.h"
#include "tangentPlane.h"

#define TINYNUM 0.00000001
#define fABS(x) ((x<0)? -x : x)
#define fEQ(x,y) (fABS(x-y)<TINYNUM)

/* parseFile takes a filename as input, tests it as an obj, and grabs all of its vertices
 * It returns of Point vector of these vertices
 * An error will be thrown if file cannot be read */
std::vector<V3> parseFile(const char* filename){
  std::vector<V3> points;
  std::ifstream fin;
  fin.open(filename);
  if(!fin.good()){
    printf("Couldn't open input file!\n");
    exit(0);
  }else{
    char line[MAX_LINE_SIZE],line_type[MAX_LINE_SIZE];
    while(fin.getline(line,MAX_LINE_SIZE)){
      float x,y,z;
      sscanf(line,"%s %f %f %f",line_type,&x,&y,&z);
      printf("line of type %s is (%.3f,%.3f,%.3f)\n",line_type,x,y,z);
      if(line_type[0]=='v'){ //vertex
        V3 p;
        p.x=x;p.y=y;p.z=z;
        points.push_back(p);
      }else printf("couldn't recognize type %s\n",line_type);
    }
  }
  fin.close();
  return points;
}

void saveMesh(std::vector<V3> points, std::vector<E> edges, char* out_filename){
  std::ofstream fout;
  fout.open(out_filename);
  if(!fout.good()){
    printf("couldn't write to %s\n",out_filename);
    exit(0);
  }else{
    for(int i=0;i<points.size();i++){
      float x,y,z;
      x=points[i].x;
      y=points[i].y;
      z=points[i].z;
      fout << "v " << x << " " << y << " " << z << "\n";
    }
    
    for(int i=0;i<edges.size();i++){
      int v1,v2;
      v1=edges[i].v1;
      v2=edges[i].v2;
      fout << "e " << v1 << " " << v2 << "\n";
    }
  }
  fout.close();
}

inline void printPoint(V3 p){
  printf("(%.3f,%.3f,%.3f)\n",p.x,p.y,p.z);
}

int main(){
  //for now we can read/write from obj file to keep things simple
  //until we decide how we're actually gonna grab kinect stuff

  //step 1: parse input file/stream
  const char* in_filename = INPUT_FILE;
  std::vector<V3> vertices = parseFile(in_filename);
  int numPoints = vertices.size();
  V3* points = (V3*) malloc(sizeof(V3)*numPoints);
  bbox system;
  for(int i=0;i<numPoints;i++){
    points[i] = vertices[i];
    printPoint(points[i]);
    system.expand(points[i]);
  }
  system.print();
  //step 2: create mesh in parallel
  //step 2a: get neighborhood for each point
  //step 2b: get centroid & PCA normals based off of neighborhoods
  Plane* planes = computeTangentPlanes(points,numPoints,0.6,0.5);
  /*for(int i=0;i<numPoints;i++){
    printf("Plane %d: \n\t",i);
    printPoint(planes[i].center);
    printf("\t");printPoint(planes[i].normal);
  }*/
  //step 2c: Propogate normal directions for surface consistency
  //step 2d: Create cubes around these centers

  //step 2e: Approximate mesh based on differences between cubes
  V3 universeSize = system.max-system.min;
  float sideLength = 0.5; //cube side length
  //modify system size to fit int num of cubes in each dir
  int widthDif = universeSize.x-std::ceil(universeSize.x/sideLength)*sideLength,
     heightDif = universeSize.y-std::ceil(universeSize.y/sideLength)*sideLength,
      depthDif = universeSize.z-std::ceil(universeSize.z/sideLength)*sideLength;
  V3 addon = V3(widthDif,heightDif,depthDif);
  addon.scale(0.5); // add evenly to min and max
  system.max += addon;
  system.min -= addon;
  V3 numCubes = system.max-system.min;
  numCubes.scale(1.0/sideLength);
  if(numCubes.x==0)numCubes.x++;
  if(numCubes.y==0)numCubes.y++;
  if(numCubes.z==0)numCubes.z++;
  printf("Num cubes in each dir is:\n");
  printPoint(numCubes);

  bbox*** cubes = (bbox***) malloc(numCubes.x*sizeof(bbox**));
  std::vector<V3> newvertices;
  std::vector<Edge> edges;
  for(int i=0;i<numCubes.x;i++){
    cubes[i] = (bbox**) malloc(numCubes.y*sizeof(bbox*));
    for(int j=0;j<numCubes.y;j++){
      cubes[i][j] = (bbox*) malloc(numCubes.z*sizeof(bbox));
      for(int k=0;k<numCubes.z;k++){
        bbox cube = bbox(system.min+V3(i*sideLength,j*sideLength,k*sideLength),sideLength,sideLength,sideLength);
        cubes[i][j][k] = cube;
        printf("cube %d,%d,%d is:\n",i,j,k);cubes[i][j][k].print();
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
          else if((v[0]<=0)==(v[1]<=0)==(v[2]<=0) ||
                  (v[0]<=0)==(v[1]<=0)==(v[3]<=0) ||
                  (v[0]<=0)==(v[2]<=0)==(v[3]<=0) ||
                  (v[1]<=0)==(v[2]<=0)==(v[3]<=0)){//three match
            for(int idx=0;idx<4;idx++){
              float vm = v[idx];V3 pm = p[idx];
              float vp = v[(idx-1)%4];V3 pp = p[(idx-1)%4];
              float vn = v[(idx+1)%4];V3 pn = p[(idx+1)%4];
              if((vm<=0)==(vp<=0) || (vm<=0)==(vn<=0)) continue; //only consider odd man out
              float frac1 = fABS(vm)/fABS(vm+vp),
                    frac2 = fABS(vm)/fABS(vm+vn);
              V3 newP1 = pm+(pp*frac1);
              V3 newP2 = pm+(pn*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(Edge(newP1,newP2));
            }
          }else{//two match
            if((v[0]<=0)==(v[1]<=0)){ //horiz line
              float frac1 = fABS(v[0])/fABS(v[0]+v[3]),
                    frac2 = fABS(v[1])/fABS(v[1]+v[2]);
              V3 newP1 = p[0]+(p[3]*frac1);
              V3 newP2 = p[1]+(p[2]*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(Edge(newP1,newP2));
            }else if((v[0]<=0)==(v[3]<=0)){ //vert line
              float frac1 = fABS(v[0])/fABS(v[0]+v[1]),
                    frac2 = fABS(v[3])/fABS(v[3]+v[2]);
              V3 newP1 = p[0]+(p[1]*frac1);
              V3 newP2 = p[3]+(p[2]*frac2);
              newvertices.push_back(newP1);
              newvertices.push_back(newP2);
              edges.push_back(Edge(newP1,newP2));
            }else{ //double diagonal
              float frac1 = fABS(v[0])/fABS(v[0]+v[1]),
                    frac2 = fABS(v[0])/fABS(v[0]+v[3]),
                    frac3 = fABS(v[2])/fABS(v[2]+v[1]),
                    frac4 = fABS(v[2])/fABS(v[2]+v[3]);
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
                edges.push_back(Edge(newP1,newP3));
                edges.push_back(Edge(newP2,newP4));
              }else{
                edges.push_back(Edge(newP1,newP2));
                edges.push_back(Edge(newP3,newP4));
              }
            }
          }
        }
      }
    }
  }

  //get rid of duplicate vertices
  int numNewVerts = newvertices.size();
  
  std::vector<int> queue;
  for(int i=0;i<numNewVerts;i++){
    V3 v1 = newvertices[i];
    for(int j=i+1;j<numNewVerts;j++){
      V3 v2 = newvertices[j];
      if(!std::isnan(v1.x) && fEQ(v1.x,v2.x) && fEQ(v1.y,v2.y) && fEQ(v1.z,v2.z)){queue.push_back(j);newvertices[j] = V3(NAN,NAN,NAN);}
    }
  }
  std::sort(queue.begin(),queue.end());
  int offset=0;
  while(!queue.empty()){
    newvertices.erase(newvertices.begin()+queue[0]+(offset++));
    queue.erase(queue.begin());
  }

  //go from edges as two vertices to edges as two indices to vertices
  std::vector<E> compressed_edges;
  for(int i=0;i<edges.size();i++){
    V3 v1=edges[i].v1,
       v2=edges[i].v2;
    int v1x=-1,
        v2x=-1;
    for(int j=0;j<newvertices.size();j++){
      if(v1x<0 && fEQ(newvertices[j].x,v1.x) && fEQ(newvertices[j].y,v1.y) && fEQ(newvertices[j].z,v1.z)) v1x=j;
      if(v2x<0 && fEQ(newvertices[j].x,v2.x) && fEQ(newvertices[j].y,v2.y) && fEQ(newvertices[j].z,v2.z)) v2x=j;
      if(v1x>=0 && v2x>=0) break;
    }
    compressed_edges.push_back(E(v1x,v2x));
  }

  //for(int i=0;i<newvertices.size();i++) printPoint(newvertices[i]);
  //for(int i=0;i<compressed_edges.size();i++) printf("edge (%d,%d),",compressed_edges[i].v1,compressed_edges[i].v2);

  //step 2f: Refine mesh

  //step 3: write mesh to output file
  int slen = 0;
  while(in_filename[slen++]!='\0');
  char out_filename[slen+5];
  for(int i=0;i<3;i++) out_filename[i] = "out"[i];
  for(int i=2;i<slen-5;i++) out_filename[i+1]=in_filename[i];
  for(int i=0;i<9;i++) out_filename[slen-4+i]="_out.obj"[i]; //9 includes null terminator
  printf("saving to %s\n",out_filename);
  saveMesh(newvertices,compressed_edges,out_filename);
  return 0;
}
