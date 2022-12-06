#ifndef REGULAR_GRID_HPP
#define REGULAR_GRID_HPP

#include "common.hpp"
#include "geom.hpp"
#include "regular_grid_boundary.hpp"
#include "agrid.hpp"

class ARegularGrid: public AGrid{
public:
	// ===== Overriden from AGrid
	int n_points() const override;
	Point point(int ipoint) const override;

	// ===== Native regular grid methods
	// coordinate-wise point number
	int nx() const;
	int ny() const;
	int nz() const;
	
	// point by ijk index
	Point point_ijk(int ix, int iy=0, int iz=0) const;

	// coordinates in direction
	const std::vector<double>& xcoo() const { return _xcoo; }
	const std::vector<double>& ycoo() const { return _ycoo; }
	const std::vector<double>& zcoo() const { return _zcoo; }

	// point index converter
	int point_ijk_to_glob(int ix, int iy, int iz) const;
	std::array<int, 3> point_glob_to_ijk(int ipoint) const;

	// face index converter
	virtual int face_ijk_to_glob(int ix, int iy, int iz, int idir) const = 0;
	virtual std::array<int, 4> face_glob_to_ijk(int face) const = 0;

	// cell index converter
	int cell_ijk_to_glob(int ix, int iy, int iz) const;
	std::array<int, 3> cell_glob_to_ijk(int icell) const;

	// steps
	double hx(int ix) const;
	double hy(int iy) const;
	double hz(int iz) const;

	// ==== boundaries
	// returns boundary for the specified type
	const RegularGridBoundary& reg_boundary(int ibnd) const;

	// defines boundary
	void define_boundary(int btype, DirectionCode dircode);
	void define_boundary(int btype, DirectionCode dircode, std::function<bool(Point)> filter);
	void define_boundary(int btype, std::shared_ptr<AGridBoundary> boundary) override;

	static std::shared_ptr<ARegularGrid> build(
		int nx, double len_x,
		int ny=1, double len_y=0,
		int nz=1, double len_z=0);
protected:
	// build grid with variable steps
	ARegularGrid(const std::vector<double>& xcoo,
	             const std::vector<double>& ycoo,
	             const std::vector<double>& zcoo);
private:
	const std::vector<double> _xcoo;
	const std::vector<double> _ycoo;
	const std::vector<double> _zcoo;
};

class RegularGrid1: public ARegularGrid{
public:
	RegularGrid1(double len_x, int nx);
	RegularGrid1(const std::vector<double>& xcoo);
	int n_faces() const override;
	int n_cells() const override;
	int n_boundary_faces() const override;
	double face_area(int iface) const override;
	double cell_volume(int icell) const override;
	Vector face_normal(int iface) const override;
	std::vector<int> tab_face_point(int iface) const override;
	std::array<int, 2> tab_face_cell(int iface) const override;
	int face_ijk_to_glob(int ix, int iy, int iz, int idir) const override;
	std::array<int, 4> face_glob_to_ijk(int face) const override;
};

class RegularGrid2: public ARegularGrid{
public:
	RegularGrid2(double len_x, int nx, double len_y, int ny);
	RegularGrid2(const std::vector<double>& xcoo,
	             const std::vector<double>& ycoo);
	int n_faces() const override;
	int n_cells() const override;
	int n_boundary_faces() const override;
	double face_area(int iface) const override;
	double cell_volume(int icell) const override;
	Vector face_normal(int iface) const override;
	std::vector<int> tab_face_point(int iface) const override;
	std::array<int, 2> tab_face_cell(int iface) const override;
	int face_ijk_to_glob(int ix, int iy, int iz, int idir) const override;
	std::array<int, 4> face_glob_to_ijk(int face) const override;
};

class RegularGrid3: public ARegularGrid{
public:
	RegularGrid3(
		double len_x, int nx,
		double len_y, int ny,
		double len_z, int nz);

	RegularGrid3(const std::vector<double>& xcoo,
	             const std::vector<double>& ycoo,
	             const std::vector<double>& zcoo);
	int n_faces() const override;
	int n_cells() const override;
	int n_boundary_faces() const override;
	double face_area(int iface) const override;
	double cell_volume(int icell) const override;
	Vector face_normal(int iface) const override;
	std::vector<int> tab_face_point(int iface) const override;
	std::array<int, 2> tab_face_cell(int iface) const override;
	int face_ijk_to_glob(int ix, int iy, int iz, int idir) const override;
	std::array<int, 4> face_glob_to_ijk(int face) const override;
};


#endif
