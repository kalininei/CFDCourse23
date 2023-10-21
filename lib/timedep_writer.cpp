#include "timedep_writer.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include "ghc/filesystem.hpp"

namespace fs = ghc::filesystem;

namespace{

void create_directory(std::string path, bool purge){
	if (fs::is_directory(path)){
		if (purge){
			fs::remove_all(path);
		} else {
			return;
		}
	}
	fs::create_directory(path);
}

}


TimeDependentWriter::TimeDependentWriter(std::string stem): _stem(stem), _series_fn(stem+".vtk.series"){
	create_directory(stem, true);

	std::ofstream ofs(_series_fn);
	if (!ofs) throw std::runtime_error("Failed to open " + _series_fn + " for writing");

	ofs << "{" << std::endl;
	ofs << "  \"file-series-version\" : \"1.0\"," << std::endl;
	ofs << "  \"files\" : [" << std::endl;
	ofs << "  ]" << std::endl;
	ofs << "}" << std::endl;
}

std::string TimeDependentWriter::add(double tm){
	std::ostringstream fn;
	fn << std::setfill('0') << std::setw(8) << std::fixed << std::setprecision(4) << tm << ".vtk";
	std::string ret = _stem + '/' + fn.str();

	std::ostringstream oss;
	if (!_fileslist.empty()){
		oss << "," << std::endl;
	}
	oss << "    {\"name\": \"" << ret << "\", \"time\": " << tm << "}";
	_fileslist += oss.str();

	save_series();

	return ret;
}

void TimeDependentWriter::save_series() const{
	std::fstream ofs(_series_fn);
	if (!ofs) throw std::runtime_error("Failed to open " + _series_fn + " for writing");

	ofs << "{" << std::endl;
	ofs << "  \"file-series-version\" : \"1.0\"," << std::endl;
	ofs << "  \"files\" : [" << std::endl;
	ofs << _fileslist << std::endl;
	ofs << "  ]" << std::endl;
	ofs << "}" << std::endl;
}
