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

struct TripletSparseMatrix: public ISparseMatrix{
	TripletSparseMatrix(size_t n_rows): _n_rows(n_rows){}
	size_t n_rows() const override{
		return _n_rows;
	};
	void set_value(size_t irow, size_t icol, double value) override{
		int addr = find_by_row_col(irow, icol);
		if (addr >= 0){
			_vals[addr] = value;
		} else {
			_rows.push_back(irow);
			_cols.push_back(icol);
			_vals.push_back(value);
		}
	}
	std::vector<double> mult(
			const std::vector<double>& x) const override{
		std::vector<double> ret(n_rows(), 0);
		for (size_t irow=0; irow < n_rows(); ++irow){
				ret[irow] = mult_row(irow, x);
		}
		return ret;
	}
	double mult_row(size_t irow,
			const std::vector<double>& x) const override{
		double sum = 0;
		for (size_t icol=0; icol < n_rows(); ++icol){
				sum += value(irow, icol) * x[icol];
		}
		return sum;
	}
	std::vector<double> diagonal() const override{
		std::vector<double> ret(n_rows());
		for (size_t i=0; i<n_rows(); ++i){
			ret[i] = value(i, i);
		}
		return ret;
	}
private:
	const size_t _n_rows;
	std::vector<int> _rows;
	std::vector<int> _cols;
	std::vector<double> _vals;

	double value(size_t irow, size_t icol) const{
		int a = find_by_row_col(irow, icol);
		if (a >= 0){
			return _vals[a];
		} else {
			return 0.0;
		}
	}

	// returns >= 0 if [irow, icol] found, else -1
	int find_by_row_col(size_t irow, size_t icol) const{
		for (size_t a=0; a<_rows.size(); ++a){
			if (_rows[a] == (int)irow
			    && _cols[a] == (int)icol){
				return (int)a;
			}
		}
		return -1;
	}
};

struct ISolver{
	virtual ~ISolver() = default;
	ISolver(const ISparseMatrix* mat): _mat(mat){}

	virtual void solve(const std::vector<double>& rhs, std::vector<double>& x) = 0;
protected:
	const ISparseMatrix* _mat;
};

struct IIterativeSolver: public ISolver{
	virtual ~IIterativeSolver() = default;
	IIterativeSolver(const ISparseMatrix* mat, size_t maxit, double eps):
		ISolver(mat), _maxit(maxit), _eps(eps){}

	void set_saver(std::shared_ptr<ASpatialApproximator> appr, const std::string& fname){
		_appr = appr;
		_writer = std::make_shared<TimeDependentWriter>(fname);
	}

	void solve(const std::vector<double>& rhs, std::vector<double>& x) override{
		for (size_t it=0; it<_maxit; ++it){
			// iteration step
			iteration_step(rhs, x);

			// save solution
			if (_appr){
				std::string fname = _writer->add(it);
				_appr->vtk_save_scalars(fname, {{"data", &x}});
			}

			// residual
			std::vector<double> lhs = _mat->mult(x);
			double norm_max = 0;
			for (size_t i=0; i<_mat->n_rows(); ++i){
				double res = std::abs(lhs[i] - rhs[i]);
				norm_max = std::max(norm_max, res);
			}
			if (norm_max < _eps){
				std::cout << "converged in " << it << " iterations" << std::endl;
				break;
			}
		}
	}
private:
	const size_t _maxit;
	const double _eps;
	std::shared_ptr<ASpatialApproximator> _appr;
	std::shared_ptr<TimeDependentWriter> _writer;

	virtual void iteration_step(const std::vector<double>& rhs, std::vector<double>& x) = 0;
};

struct JacobiSolver: public IIterativeSolver{
	JacobiSolver(const ISparseMatrix* mat, size_t maxit, double eps):
		IIterativeSolver(mat, maxit, eps){}

private:
	void iteration_step(const std::vector<double>& rhs, std::vector<double>& x) override{
		std::vector<double> diag = _mat->diagonal();
		std::vector<double> x_new(x);
		for (size_t irow=0; irow<_mat->n_rows(); ++irow){
			x_new[irow] += (rhs[irow] - _mat->mult_row(irow, x))/diag[irow];
		}
		// to the next layer: x = x_new
		std::swap(x, x_new);
	}
};


void test(){
	size_t N = 10;
	double A = 1;
	double B = 2;
	double h = 1.0/(N-1);

	ISparseMatrix* mat = new TripletSparseMatrix(N);
	std::shared_ptr<ARegularGrid> grid = RegularGrid1::build(N, 1);
	std::shared_ptr<FdmApproximator> appr = FdmApproximator::build(grid);


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
	IIterativeSolver* solver = new JacobiSolver(mat, 10'000, 1e-3);
	solver->set_saver(appr, "jacobi");

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
