#include "prob/nonstationary_problem.hpp"

void ANonstationaryProblem::solve(double tend, double tstart){
	if (!_initialized){
		std::cout << "=== Initialization" << std::endl;
		_initialize();
		for (auto& m: _monitors){
			m->initialize(tstart);
		}
		_initialized=true;
	}

	std::cout << "=== Start computation at t = " << tstart << std::endl;

	_time = tstart;
	while (_time < tend - 1e-8){
		double tau = _compute_tau();
		if (_time + tau > tend){
			tau = tend - _time;
		}
		_solve_next_step(tau);
		_time += tau;

		bool forced_break = false;
		for (auto& m: _monitors){
			if (m->run(_time)) forced_break = true;
		}

		if (forced_break) break;
	}

	for (auto& m: _monitors){
		m->finalize(_time);
	}

	std::cout << "=== Stop computation at t = " << _time << std::endl;
}

double ANonstationaryProblem::last_time() const{
	return _time;
}

void ANonstationaryProblem::add_monitor(std::shared_ptr<IMonitor> m){
	_monitors.push_back(m);
}

void ANonstationaryProblem::_initialize(){
}
