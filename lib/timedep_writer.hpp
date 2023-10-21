#include <string>

class TimeDependentWriter{
public:
	TimeDependentWriter(std::string stem);
	std::string add(double tm);
private:
	const std::string _stem;
	const std::string _series_fn;
	std::string _fileslist;
	void save_series() const;
};
