std::vector<Eigen::Vector3f> parseFilePoints(const char* filename);

void saveMesh(std::vector<Eigen::Vector3f> V, std::vector<Edge> edges, const char* out_filename);
