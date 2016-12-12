Plane getTangentPlane(Eigen::Vector3f* neighbors);
Plane* computeTangentPlanes(Eigen::Vector3f* points, int numPoints);
void insertPoints(Eigen::Vector3f* points, int numPoints, CubeData*** splitData, Eigen::Vector3f maxxyz,Eigen::Vector3f minxyz, Eigen::Vector3f diff, float cubeLength, int numCubes);
float getDist(Eigen::Vector3f p, Plane* planes, int numPlanes);
