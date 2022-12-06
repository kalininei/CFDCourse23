#include "slae/csrmat.hpp"

void CsrStencil::set_data(std::vector<int>&& addr, std::vector<int>&& cols){
	_addr = std::move(addr);
	_cols = std::move(cols);
}

int CsrStencil::n_nonzero() const{
	return (int)_cols.size();
}

int CsrStencil::n_rows() const{
	return (int)_addr.size() - 1;
}

int CsrStencil::addr_index(int irow, int jcol) const{
	for (int a=_addr[irow]; a<_addr[irow+1]; ++a){
		if (_cols[a] == jcol){
			return a;
		}
	}
	throw std::runtime_error("Failed to find csr address");
}

std::vector<int> CsrStencil::addr_diag() const{
	return std::vector<int>(_addr.begin(), _addr.end()-1);
}

void CsrStencil::add_value(int irow, int jcol, double val, std::vector<double>& mat) const{
	mat[addr_index(irow, jcol)] += val;
}

void CsrStencil::matvec(const std::vector<double>& mat, const std::vector<double>& vec, std::vector<double>& ret) const{
	for (int irow=0; irow<n_rows(); ++irow){
		ret[irow] = 0;
		for (int a=_addr[irow]; a<_addr[irow+1]; ++a){
			ret[irow] += mat[a] * vec[_cols[a]];
		}
	}
}

void CsrStencil::matvec_plus(double coeff, const std::vector<double>& mat, const std::vector<double>& vec, std::vector<double>& ret) const{
	for (int irow=0; irow<n_rows(); ++irow){
		for (int a=_addr[irow]; a<_addr[irow+1]; ++a){
			ret[irow] += coeff * mat[a] * vec[_cols[a]];
		}
	}
}

double CsrStencil::matvec_irow(int irow, const std::vector<double>& mat, const std::vector<double>& vec) const{
	double ret = 0;
	for (int a=_addr[irow]; a<_addr[irow+1]; ++a){
		ret += mat[a] * vec[_cols[a]];
	}
	return ret;
}

void CsrStencil::set_unit_diagonal(int irow, std::vector<double>& mat) const{
	mat[_addr[irow]] = 1.0;
	for (int a = _addr[irow]+1; a<_addr[irow+1]; ++a){
		mat[a] = 0.0;
	}
}

CsrStencil CsrStencil::build(const std::vector<std::set<int>>& stencil_set){
	std::vector<int> addr, cols;
	addr.push_back(0);

	for (int irow=0; irow<(int)stencil_set.size(); ++irow){
		auto& s = stencil_set[irow];
		bool has_diag = false;
		cols.push_back(irow);
		for (int jcol: s){
			if (jcol == irow){
				has_diag = true;
			} else {
				cols.push_back(jcol);
			}
		}
		addr.push_back(addr.back() + s.size());
		if (!has_diag){
			addr.back() += 1;
		}

	}

	CsrStencil ret;
	ret.set_data(std::move(addr), std::move(cols));
	return ret;
}

void CsrMatrix::matvec(const std::vector<double>& vec, std::vector<double>& ret) const{
	CsrStencil::matvec(_vals, vec, ret);
}

double CsrMatrix::matvec_irow(int irow, const std::vector<double>& vec) const{
	return CsrStencil::matvec_irow(irow, _vals, vec);
}

void CsrMatrix::matvec_plus(double coeff, const std::vector<double>& vec, std::vector<double>& ret) const{
	CsrStencil::matvec_plus(coeff, _vals, vec, ret);
}

void CsrMatrix::set_unit_diagonal(int irow){
	CsrStencil::set_unit_diagonal(irow, _vals);
}

double& CsrMatrix::value(int irow, int icol){
	int ind = addr_index(irow, icol);
	return _vals[ind];
}
