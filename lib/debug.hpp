#ifndef CFDLIB_DEBUG_HPP
#define CFDLIB_DEBUG_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"

namespace dbg{

void print(const CsrStencil& mat, const std::vector<double>& val);
void print(const CsrMatrix& mat);

}

#endif
