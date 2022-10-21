class OutputFields{
public:
	void add_scalar_data(std::string caption, const std::vector<double>* data);
	void add_vector_data(std::string caption, const std::vector<Point>* data);
};
