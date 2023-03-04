#ifndef CFDLIB_UNSTRUCTURED_GRID_HPP
#define CFDLIB_UNSTRUCTURED_GRID_HPP

#include "grid/agrid.hpp"
#include "boost/geometry.hpp"

namespace bg = boost::geometry;

class UnstructuredGrid: public AGrid{
public:
	int n_points() const override;
	int n_faces() const override;
	int n_cells() const override;
	Point point(int ipoint) const override;
	int n_boundary_faces() const override;
	double face_area(int iface) const override;
	double cell_volume(int icell) const override;
	Vector face_normal(int iface) const override;
	std::vector<int> tab_face_point(int iface) const override;
	std::vector<int> tab_cell_point(int iface) const override;
	std::array<int, 2> tab_face_cell(int iface) const override;

	static std::shared_ptr<UnstructuredGrid> read_from_vtk(std::string fn);
	void vtk_save_cells(std::string fpath) const;
	void vtk_save_faces(std::string fpath) const;

	int find_cell_index(Point p) const override;
private:
	using boost_point_t = bg::model::point<double, 3, bg::cs::cartesian>;
	using rtree_ret_t = std::pair<boost_point_t, int>;
	using rtree_t = bg::index::rtree<rtree_ret_t, bg::index::quadratic<16>>;

	UnstructuredGrid(int ndim);

	rtree_t* cell_centers_rtree() const;
	bool point_in_cell(Point p, int icell) const;
	bool point_in_cell_1d(Point p, int icell) const;
	bool point_in_cell_2d(Point p, int icell) const;
	bool point_in_cell_3d(Point p, int icell) const;

	// cell-face connection info
	struct CellFaceEntry{
		// index of the face
		int face_index;
		//  1 - if face normal is directed outside the cell
		// -1 - otherwise
		int normal_direction;
	};
	// face-cell connection info
	struct FaceCellEntry{
		// face normal direction:
		//   1d - along x-axis
		//   2d - to the right from the [point0, point1] vector
		//   3d - direction of the vector product [point0, point1]X[point1, point2] vectors
		//
		//   where point0, point1, ... - consecutive face points

		// index of the cell along the face normal direction
		// -1 if no cell (boundary face)
		int right_cell = -1;
		// index of the cell against the face normal direction
		// -1 if no cell (boundary face)
		int left_cell = -1;
	};

	std::vector<Point> _points;
	std::vector<std::vector<int>> _cells;
	std::vector<std::vector<int>> _faces;
	std::vector<int> _vtk_cell_codes;
	std::vector<FaceCellEntry> _tab_face_cell;
	std::vector<int> _internal_faces;
	std::vector<int> _boundary_faces;
	std::vector<double> _face_areas;
	std::vector<Vector> _face_normals;
	std::vector<double> _cell_volumes;

	mutable std::shared_ptr<rtree_t> _cell_centers_rtree;
};






#endif
