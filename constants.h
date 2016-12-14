#define USE_OMP
#define NUM_NEIGHBORS 5 //num neighbors for MST Propogation

struct DebugInfo{
  float fileParse,planeCreation,kNN,MST,distsFound,mesh,dedup,complete;
  int cubex,cubey,cubez;
  float redundancy;
  bbox cloud;
};

