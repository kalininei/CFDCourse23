#include "common.hpp"
#include "prog_common.hpp"
#include "slae/matrix_solver.hpp"
#include "timedep_writer.hpp"
#include "appr/fdm_approximator.hpp"
#include "grid/regular_grid.hpp"
#include <vector>
#include "tictoc.hpp"

constexpr double m_pi = 3.1415926;

struct ISparseMatrix{
	virtual ~ISparseMatrix() = default;
	virtual size_t n_rows() const = 0;
	virtual void set_value(size_t irow, size_t icol, double value) = 0;
	virtual std::vector<double> mult(const std::vector<double>& x) const = 0;
	virtual double mult_row(size_t irow, const std::vector<double>& x) const = 0;
	virtual std::vector<double> diagonal() const = 0;
};

struct DenseMatrix: public ISparseMatrix{
	DenseMatrix(size_t ncols, size_t nrows): _n_cols(ncols), _n_rows(nrows), _data(ncols*nrows, 0){
	}
	size_t n_rows() const override{
		return _n_rows;
	};

	void set_value(size_t irow, size_t icol, double value) override{
		size_t k = linear_index(irow, icol);
		_data[k] = value;
	}

	std::vector<double> mult(const std::vector<double>& x) const override{
		std::vector<double> ret(_n_rows, 0);
		for (size_t i=0; i<_n_rows; ++i){
			ret[i] += mult_row(i, x);
		}
		return ret;
	}

	double mult_row(size_t irow, const std::vector<double>& x) const override{
		double ret = 0;
		const double* it = &_data[irow * _n_rows];
		for (size_t irow=0; irow < _n_rows; ++irow){
			ret += (*it) * x[irow];
			++it;
		}
		return ret;
	}

	std::vector<double> diagonal() const override{
		std::vector<double> ret(_n_rows, 0);
		for (size_t i =0; i<_n_rows; ++i){
			size_t k = linear_index(i, i);
			ret[i] = _data[k];
		}
		return ret;
	}
private:
	const size_t _n_cols;
	const size_t _n_rows;
	std::vector<double> _data;
	// (i, j) -> k
	size_t linear_index(size_t irow, size_t icol) const{
		return irow * _n_cols + icol;
	};
};

struct MyCsrMatrix: public ISparseMatrix{
	MyCsrMatrix(size_t nrows): _addr(nrows+1, 0){

	}
	size_t n_rows() const override{
		return _addr.size() - 1;
	};

	void set_value(size_t irow, size_t icol, double value) override{
		size_t ibegin = _addr[irow];
		size_t iend = _addr[irow+1];
		auto cols_begin = _cols.begin() + ibegin;
		auto cols_end = _cols.begin() + iend;
		auto it = std::lower_bound(cols_begin, cols_end, icol);
		size_t a = it - _cols.begin();
		if (it != cols_end && *it == icol){
			_vals[a] = value;
		} else {
			for (size_t i=irow+1; i<_addr.size(); ++i) _addr[i] += 1;
			_cols.insert(_cols.begin() + a, icol);
			_vals.insert(_vals.begin() + a, value);
		}
	}

	std::vector<double> mult(const std::vector<double>& x) const override{
		std::vector<double> ret(n_rows(), 0);
		for (size_t i=0; i<n_rows(); ++i){
			ret[i] += mult_row(i, x);
		}
		return ret;
	}

	double mult_row(size_t irow, const std::vector<double>& x) const override{
		double ret = 0;
		for (size_t a = _addr[irow]; a < _addr[irow+1]; ++a){
			ret += _vals[a] * x[_cols[a]];
		}
		return ret;
	}

	std::vector<double> diagonal() const override{
		std::vector<double> ret(n_rows());
		for (size_t i=0; i<ret.size(); ++i){
			ret[i] = value(i, i);
		}
		return ret;
	}

	double value(size_t irow, size_t icol) const{
		size_t ibegin = _addr[irow];
		size_t iend = _addr[irow+1];
		auto it = std::lower_bound(_cols.begin() + ibegin, _cols.begin() + iend, icol);
		if (it != _cols.begin() + iend && *it == icol){
			size_t a = it - _cols.begin();
			return _vals[a];
		} else {
			return 0;
		}
	}
private:
	std::vector<double> _vals;
	std::vector<double> _cols;
	std::vector<double> _addr;
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

	void set_save_stride(size_t stride){
		_save_stride = stride;
	}

	void solve(const std::vector<double>& rhs, std::vector<double>& x) override{
		for (size_t it=0; it<_maxit; ++it){
			// iteration step
			iteration_step(rhs, x);

			// save solution
			if (_appr && it % _save_stride == 0){
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
			std::cout << it << " " << norm_max << std::endl;
			if (norm_max < _eps){
				std::cout << "converged in " << it << " iterations" << std::endl;
				break;
			}
		}
	}
private:
	const size_t _maxit;
	const double _eps;
	size_t _save_stride = 1;
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


struct SorSolver: public IIterativeSolver{
	SorSolver(const ISparseMatrix* mat, size_t maxit, double eps, double omega):
		IIterativeSolver(mat, maxit, eps), _omega(omega){}
private:
	void iteration_step(const std::vector<double>& rhs, std::vector<double>& x) override{
		std::vector<double> diag = _mat->diagonal();
		for (size_t irow=0; irow<_mat->n_rows(); ++irow){
			x[irow] += _omega * (rhs[irow] - _mat->mult_row(irow, x))/diag[irow];
		}
	}
	const double _omega;
};

struct SeidelSolver: public SorSolver{
	SeidelSolver(const ISparseMatrix* mat, size_t maxit, double eps):
		SorSolver(mat, maxit, eps, 1){}
};


double exact_solution(double x){
	return sin(2*m_pi*x) + 0.5*sin(10*m_pi*x);
	//return x;
}

double exact_rhs(double x){
	return 4*m_pi*m_pi*sin(2*m_pi*x) + 50*m_pi*m_pi*sin(10*m_pi*x);
	//return 0;
}

void test(){
	size_t N = 50;
	double h = 1.0/(N-1);

	//ISparseMatrix* mat = new TripletSparseMatrix(N);
	//ISparseMatrix* mat = new DenseMatrix(N, N);
	ISparseMatrix* mat = new MyCsrMatrix(N);
	std::shared_ptr<ARegularGrid> grid = RegularGrid1::build(N, 1);
	std::shared_ptr<FdmApproximator> appr = FdmApproximator::build(grid);

	Tic("matrix fill");
	// ==== Fill matrix
	// 1)
	for (size_t i=1; i<N-1; ++i){
		mat->set_value(i, i, 2.0/h/h);
		mat->set_value(i, i-1, -1.0/h/h);
		mat->set_value(i, i+1, -1.0/h/h);
	}
	// 2) u[0] = 0
	mat->set_value(0, 0, 1);
	// 3) u[N-1] = 0
	mat->set_value(N-1, N-1, 1);
	Toc("matrix fill");

	// ==== Fill rhs
	std::vector<double> f(N);
	// 1)
	for (size_t i=1; i<N-1; ++i){
		double x = i * h;
		f[i] = exact_rhs(x); 
	}
	// 2)
	f[0] = exact_solution(0);
	// 3) 
	f[N-1] = exact_solution(1);

	// ================== Solution
	std::vector<double> u(N, 0);
	//IIterativeSolver* solver = new SeidelSolver(mat, 10'000, 1e-3);
	//solver->set_saver(appr, "seidel");
	IIterativeSolver* solver = new SorSolver(mat, 10'000, 1e-3, 0.9);
	solver->set_saver(appr, "sor");
	//IIterativeSolver* solver = new JacobiSolver(mat, 10'000, 1e-3);
	//solver->set_saver(appr, "jacobi");
	solver->set_save_stride(10000000);

	Tic("matrix solver");
	solver->solve(f, u);
	Toc("matrix solver");

	//// ================== Print result
	//for (size_t i=0; i<N; ++i){
	//        std::cout << i*h << " " << u[i] << std::endl;
	//}
}

void test_csr(){
	MyCsrMatrix m(3);
	std::cout << m.value(0, 0) << std::endl;
	m.set_value(0, 0, 1);
	std::cout << m.value(0, 0) << std::endl;
	m.set_value(2, 0, 2);
	std::cout << m.value(0, 0) << std::endl;
	std::cout << m.value(2, 0) << std::endl;
}

int main(){
	try{
		//test_csr();
		test();
		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << "  " << e.what() << std::endl;
	}
}
