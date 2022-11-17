#ifndef CSRMAT_HPP
#define CSRMAT_HPP

#include "common.hpp"

// Stencil data for Csr matrix
class CsrStencil{
public:
	void set_data(std::vector<int>&& addr, std::vector<int>&& cols);

	int n_nonzero() const;
	int n_rows() const;

	const std::vector<int>& addr() const { return _addr; }
	const std::vector<int>& cols() const { return _cols; }

	// returns index of [irow, jcol] entry in values vector
	// throws if there is no such entry
	int addr_index(int irow, int jcol) const;

	// ret = mat * vec
	// mat - values vector of the sparse matrix
	void matvec(const std::vector<double>& mat, const std::vector<double>& vec, std::vector<double>& ret) const;
	// ret += coeff * mat * vec
	void matvec_plus(double coeff, const std::vector<double>& mat, const std::vector<double>& vec, std::vector<double>& ret) const;
	// mat*vec value for the specified row
	double matvec_irow(int irow, const std::vector<double>& mat, const std::vector<double>& vec) const;
	// sets diagonal value of the specified row to 1, nondiagonal to 0
	void set_unit_diagonal(int irow, std::vector<double>& mat) const;
	
	// build from vector of sets
	// i-th entry of vector -- column indices in the i-th row
	static CsrStencil build(const std::vector<std::set<int>>& stencil_set);
private:
	std::vector<int> _addr = {0};
	std::vector<int> _cols;
};


// sparse matrix in csr format.
class CsrMatrix: public CsrStencil{
public:
	CsrMatrix(const CsrStencil& stencil): CsrStencil(stencil), _vals(stencil.n_nonzero(), 0.0){}
	CsrMatrix(CsrStencil&& stencil): CsrStencil(std::move(stencil)), _vals(n_nonzero(), 0.0){}

	std::vector<double>& vals() { return _vals; }
	const std::vector<double>& vals() const { return _vals; }

	// ret = mat * vec
	void matvec(const std::vector<double>& vec, std::vector<double>& ret) const;
	// ret += coeff * mat * vec
	void matvec_plus(double coeff, const std::vector<double>& vec, std::vector<double>& ret) const;
	// mat*vec value for the specified row
	double matvec_irow(int irow, const std::vector<double>& vec) const;
	// sets diagonal value of the specified row to 1, nondiagonal to 0
	void set_unit_diagonal(int irow);
	// returns value by row and column indices
	double& value(int irow, int icol);
private:
	std::vector<double> _vals;
};


#endif
