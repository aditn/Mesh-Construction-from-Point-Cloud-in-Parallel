std::vector<V3> parseFilePoints(const char* filename);

void saveMesh(std::vector<V3> V, std::vector<Edge> edges, const char* out_filename);

void savePlanes(Plane* planes, int numPoints,const char* out_filename);
