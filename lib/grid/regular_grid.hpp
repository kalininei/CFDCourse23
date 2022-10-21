#ifndef REGULAR_GRID_HPP
#define REGULAR_GRID_HPP

#include "common.hpp"
#include "geom.hpp"
#include "regular_grid_boundary.hpp"

class RegularGrid{
public:
	// geometric dimension
	const int dim;

	// build grid with constant steps
	RegularGrid(
		double len_x, int nx,
		double len_y=0, int ny=0,
		double len_z=0, int nz=0);

	// build grid with variable steps
	RegularGrid(const std::vector<double>& xcoo,
	            const std::vector<double>& ycoo,
	            const std::vector<double>& zcoo);

	// sizes
	int n_points() const;

	// coordinate-wise sizes
	int nx() const;
	int ny() const;
	int nz() const;

	// point by the global index
	Point point(int point_index) const;

	// point by ijk index
	Point point_ijk(int ix, int iy=0, int iz=0) const;

	const std::vector<double>& xcoo() const { return _xcoo; }
	const std::vector<double>& ycoo() const { return _ycoo; }
	const std::vector<double>& zcoo() const { return _zcoo; }

	// point index converter
	int ijk_to_glob(int ix, int iy, int iz) const;
	std::array<int, 3> glob_to_ijk(int point_index) const;

	// steps
	double hx(int ix) const;
	double hy(int iy) const;
	double hz(int iz) const;

	// ==== boundaries
	// passed boundary types
	std::vector<int> btypes() const;
	
	// returns boundary for the specified type
	const RegularGridBoundary& boundary(int ibnd) const;

	// defines boundary
	void define_boundary(int btype, DirectionCode dircode);
	void define_boundary(int btype, DirectionCode dircode, std::function<bool(Point)> filter);
private:
	const std::vector<double> _xcoo;
	const std::vector<double> _ycoo;
	const std::vector<double> _zcoo;

	std::map<int, std::shared_ptr<RegularGridBoundary>> _boundary;
};


#endif
