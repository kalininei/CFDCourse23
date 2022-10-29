#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <fstream>

constexpr double PI2 = 8*std::atan(1.0);

double exact_solution(double x){
	return sin(PI2*2*x) + 0.4*cos(PI2*7*x);
}

double exact_solution_d1(double x){
	return 2*PI2*cos(PI2*2*x) - 0.4*7*PI2*sin(PI2*7*x);
}

std::vector<double> tabulate(int n_nodes, std::function<double(double)> analytic_func){
	std::vector<double> ret(n_nodes);

	for (int i=0; i<n_nodes; ++i){
		double x = double(i)/(n_nodes-1);
		ret[i] = analytic_func(x);
	}

	return ret;
}

double integrate01(const std::vector<double>& u){
	double h = 1.0 / (u.size() - 1);
	double ret = 0;

	for (size_t i=0; i<u.size()-1; ++i){
		ret += (u[i] + u[i+1])/2 * h;
	}

	return ret;
}

double norm2(const std::vector<double>& u){
	double volume = 1.0;
	std::vector<double> u2(u.size());
	for (size_t i=0; i<u2.size(); ++i){
		u2[i] = u[i] * u[i];
	}
	double integral = integrate01(u2);

	return std::sqrt(integral) / volume;
}
double norm2(const std::vector<double>& u1, const std::vector<double>& u2){
	std::vector<double> diff(u1.size()); // u1 - u2
	for (size_t i=0; i<u1.size(); ++i){
		diff[i] = u1[i] - u2[i];
	}
	return norm2(diff);
}

std::vector<double> euler_solution(int n_nodes, double y0, std::function<double(double)> rhs){
	std::vector<double> ret(n_nodes);
	ret[0] = y0;
	double h = 1.0 / (n_nodes-1);
	for (int i=0; i<n_nodes-1; ++i){
		double x = i*h;
		ret[i+1] = ret[i] + h*rhs(x);
	}
	return ret;
}

int main(){
	try{
		for (int n_cells: {10, 20, 50, 100, 200, 500, 1000, 2000, 10'000, 20'000, 100'000}){
			int n_nodes = n_cells + 1;

			std::vector<double> exact = tabulate(n_nodes, exact_solution);
			std::vector<double> numer = euler_solution(n_nodes, exact_solution(0), exact_solution_d1);

			double euler_error = norm2(exact, numer);
			std::cout << n_cells << " " << euler_error << std::endl;
		}
	} catch (const std::runtime_error& e){
		std::cerr << "Exception caught: " << e.what() << std::endl;
	}
	std::cout << "DONE" << std::endl;
}
