#ifndef _DENSE_MAT_HPP_
#define _DENSE_MAT_HPP_

#include <vector>
#include <array>
#include "geom.hpp"

// returns full dense matrix from given lower diagonal values
// lower is given as a plain array without upper zeros
std::vector<double> dense_from_lower(const std::vector<double>& lower);

// determinant of 3x3 matrix
double determinant3(const std::array<double, 9>& mat);

// mat => |mat| * tranpose(invert(mat))
std::array<double, 9> det_transpose_invert3(const std::array<double, 9>& mat);

// => {mat3x3 * point_i}
std::vector<Point> dense_mult_point(const std::array<double, 9>& mat, const std::vector<Point>& points);

// => v[i]*v[j], only for i<=j
std::vector<double> outer_lower(const std::vector<double>& v);

// => dot(v[i], v[j]), only for i<=j
std::vector<double> vector_dot_outer_lower(const std::vector<Vector>& v);



#endif
