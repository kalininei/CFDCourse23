#include "common.hpp"
#include "prog_common.hpp"
#include "slae/matrix_solver.hpp"
#include "timedep_writer.hpp"
#include "appr/fdm_approximator.hpp"
#include "grid/regular_grid.hpp"
#include <vector>

struct ISparseMatrix{
	virtual ~ISparseMatrix() = default;
	virtual size_t n_rows() const = 0;
	virtual void set_value(size_t irow, size_t icol, double value) = 0;
	virtual std::vector<double> mult(const std::vector<double>& x) const = 0;
	virtual double mult_row(size_t irow, const std::vector<double>& x) const = 0;
	virtual std::vector<double> diagonal() const = 0;
};

struct ISolver{
	virtual ~ISolver() = default;
	ISolver(const ISparseMatrix* mat): _mat(mat){}

	virtual void solve(const std::vector<double>& rhs, std::vector<double>& x) = 0;
protected:
	const ISparseMatrix* _mat;
};

struct JacobiSolver: public ISolver{
	JacobiSolver(const ISparseMatrix* mat, size_t maxit, double eps):
		ISolver(mat), _maxit(maxit), _eps(eps){}

	void solve(const std::vector<double>& rhs, std::vector<double>& x) override{
		std::vector<double> diag = _mat->diagonal();
		for (size_t it=0; it<_maxit; ++it){
			// iterations
			std::vector<double> x_new(x);
			for (size_t irow=0; irow<_mat->n_rows(); ++irow){
				x_new[irow] += (rhs[irow] - _mat->mult_row(irow, x))/diag[irow];
			}

			// to the next layer: x = x_new
			std::swap(x, x_new);

			// residual
			std::vector<double> lhs = _mat->mult(x);
			double norm_max = 0;
			for (size_t i=0; i<_mat->n_rows(); ++i){
				double res = std::abs(lhs[i] - rhs[i]);
				norm_max = std::max(norm_max, res);
			}
			if (norm_max < _eps){
				break;
			}
		}
	}
private:
	const size_t _maxit;
	const double _eps;
};


void test(){
	ISparseMatrix* mat;
	size_t N = 10;
	double A = 1;
	double B = 2;
	double h = 1.0/(N-1);

	std::shared_ptr<ARegularGrid> grid = RegularGrid1::build(N, 1);
	std::shared_ptr<FdmApproximator> appr = FdmApproximator::build(grid);


	std::vector<double> data(N, 0);
	TimeDependentWriter writer("jacobi");
	{
		std::string fname = writer.add(0);
		appr->vtk_save_scalars(fname, {{"data", &data}});
	}
	{
		std::string fname = writer.add(1);
		appr->vtk_save_scalars(fname, {{"data", &data}});
	}

	// ==== Fill matrix
	// 1)
	for (size_t i=1; i<N-1; ++i){
		mat->set_value(i, i, 2.0/h/h);
		mat->set_value(i, i-1, -1.0/h/h);
		mat->set_value(i, i+1, -1.0/h/h);
	}
	// 2) u[0] = A
	mat->set_value(0, 0, 1);
	// 3) u[N-1] = B
	mat->set_value(N-1, N-1, 1);

	// ==== Fill rhs
	std::vector<double> f(N);
	// 1)
	for (size_t i=1; i<N-1; ++i){
		f[i] = 0;
	}
	// 2)
	f[0] = A;
	// 3) 
	f[N-1] = B;

	// ================== Solution
	std::vector<double> u(N, 0);
	ISolver* solver = new JacobiSolver(mat, 10'000, 1e-3);

	solver->solve(f, u);

	// ================== Print result
	for (size_t i=0; i<N; ++i){
		std::cout << i*h << " " << u[i] << std::endl;
	}
}

int main(){
	try{
		test();
		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << "  " << e.what() << std::endl;
	}
}
