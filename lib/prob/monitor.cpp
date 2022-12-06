#include <sstream>
#include <fstream>
#include <iomanip>
#include <list>
#include "prob/monitor.hpp"

AMonitor_TimeGap::AMonitor_TimeGap(double delta_t, bool from_zero): IMonitor(), _delta_t(delta_t), _from_zero(from_zero){ }

void AMonitor_TimeGap::initialize(double tstart){
	_tstart = tstart;
	_run_calls = 0;
	_apply_calls = 0;
	_t0 = (_from_zero) ? 0 : tstart;
	_last_apply_gap = 0;

	_initialize_core(tstart);
}

bool AMonitor_TimeGap::run(double tcurrent){
	bool ret = false;
	tcurrent -= _t0;
	int igap = tcurrent/_delta_t;
	if (_delta_t - (tcurrent -_delta_t*igap) < 1e-8){
		igap += 1;
	}
	if (igap > _last_apply_gap){
		ret = _apply(tcurrent);
		_last_apply_gap = igap;
		_apply_calls += 1;
	}
	_run_calls += 1;
	return ret;
}

void AMonitor_TimeGap::finalize(double tend){ }

AMonitor_IterGap::AMonitor_IterGap(int delta_iter, int start_iter): IMonitor(), _delta_iter(delta_iter), _start_iter(start_iter){ }

void AMonitor_IterGap::initialize(double tstart){
	_run_calls = 0;
	_apply_calls = 0;
	_initialize_core(tstart);
}

bool AMonitor_IterGap::run(double tcurrent){
	bool ret = false;
	int cur_iter = _start_iter + _run_calls + 1;

	if (cur_iter % _delta_iter == 0){
		ret = _apply(tcurrent, cur_iter);
		_apply_calls += 1;
	}
	_run_calls += 1;
	return ret;
}

void AMonitor_IterGap::finalize(double tend){ }

ConsoleIterReport::ConsoleIterReport(int delta_iter, int start_iter): AMonitor_IterGap(delta_iter, start_iter){ }

void ConsoleIterReport::_initialize_core(double tstart){ }

bool ConsoleIterReport::_apply(double tcurrent, int cur_iter){
	std::ostringstream oss;

	oss << "step=" << std::setw(5) << cur_iter << ", ";
	oss << "t=" << std::fixed << std::setprecision(3) << tcurrent << " / " << _additional_info();

	std::cout << oss.str() << std::endl;
	return false;
}

std::string ConsoleIterReport::_additional_info(){
	return "";
}

AVtkFieldSaver::AVtkFieldSaver(std::string filename_stem, const ASpatialApproximator* appr):
		_filename_stem(filename_stem),
		_appr(appr){
	size_t pos = filename_stem.rfind("/");
	if (pos == filename_stem.npos){
		_filename_startpos = 0;
	} else {
		_filename_startpos = (int)pos+1;
	}
}

void AVtkFieldSaver::add_fun(std::string dataname, func_t func){
	_data.emplace_back(dataname, func);
}

void AVtkFieldSaver::_initialize_saver(double tstart){
	_written.clear();
	_last_written_iter = -1;
}

void AVtkFieldSaver::_finalize_saver(double tend, int cur_iter){
	if (cur_iter != _last_written_iter){
		_apply_saver(tend, cur_iter);
	}
}

void AVtkFieldSaver::_write_series() const{
	std::string fn = _filename_stem + ".vtk.series";
	std::ofstream ofs(fn);
	if (!ofs) throw std::runtime_error("Failed to open " + fn + " for writing");

	ofs << "{" << std::endl;
	ofs << "  \"file-series-version\" : \"1.0\"," << std::endl;
	ofs << "  \"files\" : [" << std::endl;
	for (size_t i=0; i<_written.size(); ++i){
		const auto& it = _written[i];
		ofs << "    {\"name\": \"" << it.first << "\", \"time\": " << it.second << "}";
		if (i!=_written.size()-1){
			ofs << ",";
		}
		ofs << std::endl;
	}
	ofs << "  ]" << std::endl;
	ofs << "} " << std::endl;
}

bool AVtkFieldSaver::_apply_saver(double tcurrent, int cur_iter){
	std::ostringstream fn;
	fn << _filename_stem << "_" << std::setfill('0') << std::setw(6) << cur_iter << ".vtk";
	std::cout << "Save data into " << fn.str() << std::endl;

	std::list<std::vector<double>> data;
	std::map<std::string, const std::vector<double>*> scalars;

	for (auto& it: _data){
		data.push_back(it.second());
		scalars[it.first] = &data.back();
	}

	_appr->vtk_save_scalars(fn.str(), scalars);

	_written.emplace_back(fn.str().substr(_filename_startpos), tcurrent);
	_last_written_iter = cur_iter;
	_write_series();
	return false;
}

VtkFieldIterSaver::VtkFieldIterSaver(
		int delta_iter,
		std::string filename_stem, 
		const ASpatialApproximator* appr,
		int start_iter): AMonitor_IterGap(delta_iter, start_iter),
		                 AVtkFieldSaver(filename_stem, appr){
}


void VtkFieldIterSaver::_initialize_core(double tstart){
	_initialize_saver(tstart);
}

bool VtkFieldIterSaver::_apply(double tcurrent, int cur_iter){
	return _apply_saver(tcurrent, cur_iter);
}

void VtkFieldIterSaver::finalize(double tend){
	_finalize_saver(tend, _start_iter + _run_calls);
}

VtkFieldTimeSaver::VtkFieldTimeSaver(
		double delta_t,
		std::string filename_stem, 
		const ASpatialApproximator* appr,
		bool from_zero): AMonitor_TimeGap(delta_t, from_zero),
		                 AVtkFieldSaver(filename_stem, appr){
}


void VtkFieldTimeSaver::_initialize_core(double tstart){
	_initialize_saver(tstart);
}

bool VtkFieldTimeSaver::_apply(double tcurrent){
	return _apply_saver(tcurrent, _run_calls+1);
}

void VtkFieldTimeSaver::finalize(double tend){
	_finalize_saver(tend, _run_calls);
}

FunctionalSaver::FunctionalSaver(
		double delta_iter,
		std::string filename,
		int start_iter): AMonitor_IterGap(delta_iter, start_iter), _filename(filename){ }

void FunctionalSaver::add_fun(std::string caption, func_t fun){
	_funcs[caption] = fun;
}

void FunctionalSaver::_initialize_core(double tstart){
	_ofs.open(_filename);
	_ofs << std::setw(16) << "Time step "
		<< std::setw(16) << "Time ";
	for (auto& it: _funcs){
		_ofs << std::setw(16) << it.first << " ";
	}
	_ofs << std::endl;
}

bool FunctionalSaver::_apply(double tcurrent, int cur_iter){
	_ofs << std::setw(16) << cur_iter << " "
		<< std::setw(16) << std::scientific << tcurrent;
	for (auto& it: _funcs){
		_ofs << std::setw(16) << std::scientific << it.second() << " ";
	}
	_ofs << std::endl;
	return false;
}

void FunctionalSaver::finalize(double tend){
	if (_ofs.is_open()) _ofs.close();
}
