#ifndef MONITOR_HPP
#define MONITOR_HPP

#include "common.hpp"
#include <fstream>
#include "appr/spatial_approximator.hpp"

class IMonitor{
public:
	virtual ~IMonitor() = default;

	virtual void initialize(double tstart) = 0;
	virtual bool run(double tcurrent) = 0;
	virtual void finalize(double tend) = 0;
};

class AMonitor_TimeGap: public IMonitor{
public:
	AMonitor_TimeGap(double delta_t, bool from_zero=true);
	virtual ~AMonitor_TimeGap() = default;

	void initialize(double tstart) override;
	bool run(double tcurrent) override;
	void finalize(double tend) override;
protected:
	const double _delta_t;
	const bool _from_zero;

	double _t0 = 0;
	double _tstart = 0;
	int _last_apply_gap = 0;
	int _run_calls = 0;
	int _apply_calls = 0;

	virtual void _initialize_core(double tstart) = 0;
	virtual bool _apply(double tcurrent) = 0;
};

class AMonitor_IterGap: public IMonitor{
public:
	AMonitor_IterGap(int delta_iter, int start_iter=0);
	virtual ~AMonitor_IterGap() = default;

	void initialize(double tstart) override;
	bool run(double tcurrent) override;
	void finalize(double tend) override;
protected:
	const int _delta_iter;
	const int _start_iter;
	int _run_calls = 0;
	int _apply_calls = 0;

	virtual void _initialize_core(double tstart) = 0;
	virtual bool _apply(double tcurrent, int cur_iter) = 0;
};

class ConsoleIterReport: public AMonitor_IterGap{
public:
	ConsoleIterReport(int delta_iter, int start_iter=0);
protected:
	void _initialize_core(double tstart) override;
	bool _apply(double tcurrent, int cur_iter) override;
	virtual std::string _additional_info();
};

class AVtkFieldSaver{
public:
	using func_t = std::function<std::vector<double>()>;

	AVtkFieldSaver(std::string filename_stem, const ASpatialApproximator* appr);
	virtual ~AVtkFieldSaver() = default;

	void add_fun(std::string dataname, func_t func);
protected:
	void _initialize_saver(double tstart);
	bool _apply_saver(double tcurrent, int cur_iter);
	void _finalize_saver(double tend, int cur_iter);
private:
	std::vector<std::pair<std::string, func_t>> _data;
	const std::string _filename_stem;
	const ASpatialApproximator* _appr;
	int _filename_startpos;
	std::vector<std::pair<std::string, double>> _written;
	int _last_written_iter = 0;

	void _write_series() const;
};

class VtkFieldIterSaver: public AMonitor_IterGap, public AVtkFieldSaver{
public:
	VtkFieldIterSaver(
		int delta_iter,
		std::string filename_stem, 
		const ASpatialApproximator* appr,
		int start_iter=0);

	void finalize(double tend) override;
protected:
	void _initialize_core(double tstart) override;
	bool _apply(double tcurrent, int cur_iter) override;
};

class VtkFieldTimeSaver: public AMonitor_TimeGap, public AVtkFieldSaver{
public:
	VtkFieldTimeSaver(
		double delta_time,
		std::string filename_stem, 
		const ASpatialApproximator* appr,
		bool from_zero=true);

	void finalize(double tend) override;
protected:
	void _initialize_core(double tstart) override;
	bool _apply(double tcurrent) override;
};

class FunctionalSaver: public AMonitor_IterGap{
public:
	using func_t = std::function<double()>;

	FunctionalSaver(
		double delta_iter,
		std::string filename,
		int start_iter=0);

	void add_fun(std::string caption, func_t);
protected:
	std::map<std::string, func_t> _funcs;
	std::ofstream _ofs;
	std::string _filename;

	void _initialize_core(double tstart) override;
	bool _apply(double tcurrent, int cur_iter) override;
	void finalize(double tend) override;
};

#endif
