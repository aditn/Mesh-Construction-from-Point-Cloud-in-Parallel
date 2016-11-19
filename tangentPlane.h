Plane getTangentPlane(std::vector<V3> neighbors);
Plane* computeTangentPlanes(V3* points, int numPoints, float ro, float delta);
float getDist(V3 p, Plane* planes, int numPlanes);
