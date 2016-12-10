std::vector<Eigen::Vector3f> parseFilePoints(const char* filename);

void savePlanes(Plane* planes, int numPoints,const char* out_filename);
void saveMesh(std::vector<Eigen::Vector3f> V, std::vector<Edge> edges, const char* out_filename);
