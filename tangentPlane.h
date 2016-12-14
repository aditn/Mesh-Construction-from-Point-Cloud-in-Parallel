Plane* computeTangentPlanes(Eigen::Vector3f* points, int numPoints);

std::vector<Edge> getNeighborEdges(Eigen::Vector3f* point,int numPoints,Plane* planes);  
float getDist(Eigen::Vector3f p, Plane* planes, int numPlanes);
float getDistLookup(Eigen::Vector3f p,Plane* planes,int numPlanes);
void setDists(float* dists, Eigen::Vector3f* ps,int numPoints, Plane* planes, int numPlanes);
