Plane* computeTangentPlanes(Eigen::Vector3f* points, int numPoints);

std::vector<int> getNearest(Eigen::Vector3f* points,int numPoints, int idx, int numNeighbors);

float getDist(Eigen::Vector3f p, Plane* planes, int numPlanes);

void setDists(float* dists, Eigen::Vector3f* ps,int numPoints, Plane* planes, int numPlanes);
