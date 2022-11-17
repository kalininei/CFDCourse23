#include "common.hpp"
#include "prog_common.hpp"
#include "slae/matrix_solver.hpp"

constexpr double PI2 = 8*std::atan(1.0);

double exact_solution(double x){
	return sin(PI2*2*x) + 0.4*cos(PI2*7*x);
}

double exact_solution_d1(double x){
	return 2*PI2*cos(PI2*2*x) - 0.4*7*PI2*sin(PI2*7*x);
}

double exact_solution_d2(double x){
	return -2*PI2*2*PI2*sin(PI2*2*x) - 0.4*7*7*PI2*PI2*cos(PI2*7*x);
}

double rhs(double x){
	return -exact_solution_d2(x);
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

std::vector<double> fdm_solution_poisson(int n_nodes, double left_value, double right_value,
                                         std::function<double(double)> rhs_fun){
	double h = 1.0 / (n_nodes - 1);
	// fill stencil
	std::vector<std::set<int>> stencil_set(n_nodes);
	for (int irow=0; irow < n_nodes; ++irow){
		std::set<int> cols;
	
		// main diagonal
		cols.insert(irow);
		// lower diagonal
		if (irow > 0) cols.insert(irow - 1);
		// upper diagonal
		if (irow < n_nodes - 1) cols.insert(irow + 1);

		stencil_set[irow] = cols;
	}
	CsrStencil stencil = CsrStencil::build(stencil_set);

	// fill values
	CsrMatrix matrix(stencil);

	// boundary conditions
	matrix.value(0, 0) = 1;
	matrix.value(n_nodes-1, n_nodes-1) = 1;

	// internal
	for (int irow=1; irow < n_nodes-1; ++irow){
		matrix.value(irow, irow-1) = -1.0/h/h;
		matrix.value(irow, irow+1) = -1.0/h/h;
		matrix.value(irow, irow) = 2.0/h/h;
	}

	// rhs assembling
	std::vector<double> rhs_tab = tabulate(n_nodes, rhs_fun);
	rhs_tab[0] = left_value;
	rhs_tab[n_nodes-1] = right_value;

	// solve matrix
	AmgcMatrixSolver slv;
	slv.set_matrix(matrix);

	std::vector<double> ret(n_nodes, 0);
	slv.solve(rhs_tab, ret);

	return ret;
}


int main(){
	try{
		for (int n_cells: {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10'000}){
			int n_nodes = n_cells + 1;
			std::vector<double> exact = tabulate(n_nodes, exact_solution);
			std::vector<double> numer = fdm_solution_poisson(n_nodes, exact_solution(0), exact_solution(1), rhs);

			double error = norm2(exact, numer);
			std::cout << n_cells << " " << error << std::endl;
		}

		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << "  " << e.what() << std::endl;
	}
}
