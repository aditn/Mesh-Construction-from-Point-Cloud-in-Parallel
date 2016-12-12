Plane getTangentPlane(Eigen::Vector3f* neighbors);
Plane* computeTangentPlanes(Eigen::Vector3f* points, int numPoints);
void insertPoints(Eigen::Vector3f* points, int numPoints, CubeData*** splitData, Eigen::Vector3f maxxyz,Eigen::Vector3f minxyz, Eigen::Vector3f diff, float cubeLength, int numCubes);
std::vector<int> getNearest(Eigen::Vector3f* points,int numPoints, int idx, int numNeighbors);
float getDist(Eigen::Vector3f p, Plane* planes, int numPlanes);
void setDists(float* dists, Eigen::Vector3f* ps,int numPoints, Plane* planes, int numPlanes);
