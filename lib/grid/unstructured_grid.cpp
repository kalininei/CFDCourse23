#include "unstructured_grid.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>

namespace{

// valid cell types
enum struct CellCode{
	SEGMENT = 3,
	POLYGON = 7,
	TETRAHEDRON = 10,
	HEXAHEDRON = 12,
	WEDGE = 13,
	PYRAMID = 14,
	PENTAPRISM = 15,
	HEXAPRISM = 16
};


struct AVtkCell{
	AVtkCell(std::vector<int>&& p): points(std::move(p)) {}
	virtual ~AVtkCell() = default;

	// all faces normals are outward with respect to the cell
	std::vector<std::vector<int>> faces() const{
		std::vector<std::vector<int>> ret = faces_local();
		for (auto& it1: ret)
		for (auto& it2: it1){
			it2 = points[it2];
		}
		return ret;
	}
	virtual std::vector<std::vector<int>> faces_local() const = 0;

	virtual int dim() const = 0;
	virtual int code() const = 0;

	std::vector<int> points;
	static std::shared_ptr<AVtkCell> build(std::vector<int> points, int code);
};

struct VertexCell: public AVtkCell{
	VertexCell(std::vector<int>&& p): AVtkCell(std::move(p)) {}

	std::vector<std::vector<int>> faces_local() const override{
		throw std::runtime_error("Vertex cells has no faces");
	}

	int dim() const override { return 0; }
	int code() const override { return 1; }
};

struct LineCell: public AVtkCell{
	LineCell(std::vector<int>&& p): AVtkCell(std::move(p)) {}
	std::vector<std::vector<int>> faces_local() const override{
		return { {0}, {1} };
	}
	int dim() const override { return 1; }
	int code() const override { return 3; }
};

struct PolygonCell: public AVtkCell{
	PolygonCell(std::vector<int>&& p): AVtkCell(std::move(p)) {}
	std::vector<std::vector<int>> faces_local() const override{
		std::vector<std::vector<int>> ret;
		for (int i=0; i<(int)points.size(); ++i){
			int inext = (i==(int)points.size()-1) ? 0 : i+1;
			ret.push_back({i, inext});
		}
		return ret;
	}
	int dim() const override { return 2; }
	int code() const override { return 7; }
};

struct TetraCell: public AVtkCell{
	TetraCell(std::vector<int>&& p): AVtkCell(std::move(p)) {}
	std::vector<std::vector<int>> faces_local() const override{
		return {
			{0, 2, 1},
			{1, 2, 3},
			{0, 1, 3},
			{0, 3, 2}
		};
	}
	int dim() const override { return 3; }
	int code() const override { return 10; }
};

struct PyramidCell: public AVtkCell{
	PyramidCell(std::vector<int>&& p): AVtkCell(std::move(p)) {}
	std::vector<std::vector<int>> faces_local() const override{
		return {
			{0, 3, 2, 1},
			{0, 4, 3},
			{2, 3, 4},
			{1, 4, 0},
			{1, 2, 4}
		};
	}
	int dim() const override { return 3; }
	int code() const override { return 14; }
};

struct WedgeCell: public AVtkCell{
	WedgeCell(std::vector<int>&& p): AVtkCell(std::move(p)) {}
	std::vector<std::vector<int>> faces_local() const override{
		return {
			{0, 3, 4, 1},
			{3, 5, 4},
			{0, 1, 2},
			{1, 4, 5, 2},
			{0, 2, 5, 3}
		};
	}
	int dim() const override { return 3; }
	int code() const override { return 13; }
};

struct HexaCell: public AVtkCell{
	HexaCell(std::vector<int>&& p): AVtkCell(std::move(p)) {}
	std::vector<std::vector<int>> faces_local() const override{
		return {
			{0, 1, 5, 4},
			{3, 7, 6, 2},
			{1, 2, 6, 5},
			{0, 4, 7, 3},
			{0, 3, 2, 1},
			{4, 5, 6, 7}
		};
	}
	int dim() const override { return 3; }
	int code() const override { return 12; }
};

struct PentaPrismCell: public AVtkCell{
	PentaPrismCell(std::vector<int>&& p): AVtkCell(std::move(p)) {}
	std::vector<std::vector<int>> faces_local() const override{
		return {
			{0, 4, 3, 2, 1},
			{5, 6, 7, 8, 9},
			{0, 1, 6, 5},
			{1, 2, 7, 6},
			{2, 3, 8, 7},
			{3, 4, 9, 8},
			{4, 0, 5, 9}
		};
	}
	int dim() const override { return 3; }
	int code() const override { return 15; }
};

struct HexaPrismCell: public AVtkCell{
	HexaPrismCell(std::vector<int>&& p): AVtkCell(std::move(p)) {}

	std::vector<std::vector<int>> faces_local() const override{
		return {
			{0, 5, 4, 3, 2, 1},
			{6, 7, 8, 9, 10, 11},
			{0, 1, 7, 6},
			{1, 2, 8, 7},
			{2, 3, 9, 8},
			{3, 4, 10, 9},
			{4, 5, 11, 10},
			{5, 0, 6, 11}
		};
	}

	int dim() const override { return 3; }
	int code() const override { return 16; }
};

std::shared_ptr<AVtkCell> AVtkCell::build(std::vector<int> points, int code){
	if (code == 1){
		return std::make_shared<VertexCell>(std::move(points));
	} else if (code == 3){
		return std::make_shared<LineCell>(std::move(points));
	} else if (code == 5 || code == 7 || code == 9){
		return std::make_shared<PolygonCell>(std::move(points));
	} else if (code == 8){
		std::swap(points[2], points[3]);
		return std::make_shared<PolygonCell>(std::move(points));
	} else if (code == 10){
		return std::make_shared<TetraCell>(std::move(points));
	} else if (code == 12){
		return std::make_shared<HexaCell>(std::move(points));
	} else if (code == 11){
		std::swap(points[2], points[3]);
		std::swap(points[6], points[7]);
		return std::make_shared<HexaCell>(std::move(points));
	} else if (code == 13){
		return std::make_shared<WedgeCell>(std::move(points));
	} else if (code == 14){
		return std::make_shared<PyramidCell>(std::move(points));
	} else if (code == 15){
		return std::make_shared<PentaPrismCell>(std::move(points));
	} else if (code == 16){
		return std::make_shared<HexaPrismCell>(std::move(points));
	} else {
		throw std::runtime_error("Unsupported vtk cell type: " + std::to_string(code));
	}
}


struct VtkFace{
	static constexpr int MAX_FACE_LEN = 6;

	VtkFace(const std::vector<int>& vpoints, bool& was_reverted){
		if (vpoints.size() > MAX_FACE_LEN){
			throw std::runtime_error("Faces with more than " 
			                         + std::to_string(MAX_FACE_LEN) 
			                         + " are not supported");
		}
		was_reverted = false;
		points_size = vpoints.size();
		std::copy(vpoints.begin(), vpoints.end(), points.begin());
		std::fill(points.begin()+points_size, points.end(), -1);

		if (vpoints.size() == 2){
			if (points[0] > points[1]){
				std::swap(points[0], points[1]);
				was_reverted = true;
			}
		} else if (vpoints.size() > 2){
			auto minit = std::min_element(points.begin(), points.begin() + points_size);
			std::rotate(points.begin(), minit, points.begin() + points_size);
			if (points[1] > points[points_size-1]){
				was_reverted = true;
				std::reverse(points.begin()+1, points.begin()+points_size);
			}
		}
	};

	bool is_boundary() const{
		return ((left_cell<0 && right_cell>=0)
		       || (left_cell>=0 && right_cell<0));
	};

	std::vector<int> points_vec() const{
		return std::vector<int>(points.begin(), points.begin() + points_size);
	};

	std::array<int, MAX_FACE_LEN> points;
	int index = -1;
	int points_size = 0;
	int left_cell = -1;
	int right_cell = -1;

	friend bool operator<(const VtkFace& f1, const VtkFace& f2){
		return f1.points < f2.points;
	};
};


struct VtkFaceData{
	static constexpr int MAX_FACE_LEN = 6;
	std::set<VtkFace> faces;

	int n_faces() const {
		return (int)faces.size();
	}

	std::vector<std::vector<int>> faces_to_vec() const{
		std::vector<std::vector<int>> ret;
		for (auto& it: faces){
			ret.push_back(it.points_vec());
		}
		return ret;
	}

	void add_face(const std::vector<int>& face, int icell){
		bool was_reverted;
		auto ires = faces.insert(VtkFace(face, was_reverted));
		VtkFace& fc = const_cast<VtkFace&>(*ires.first);
		if (was_reverted){
			fc.right_cell = icell;
		} else {
			fc.left_cell = icell;
		}
	}

	const VtkFace* find_face(const std::vector<int>& points){
		bool r;
		VtkFace f(points, r);
		auto fnd = faces.find(f);
		if (fnd == faces.end()){
			return nullptr;
		} else {
			return &(*fnd);
		}
	}
};

VtkFaceData assemble_faces_vtk(const std::vector<std::shared_ptr<AVtkCell>>& cells){
	VtkFaceData ret;
	for (int icell=0; icell<(int)cells.size(); ++icell){
		for (const std::vector<int>& face: cells[icell]->faces()){
			ret.add_face(face, icell);
		}
	}
	int index = 0;
	for (auto& f: ret.faces){
		VtkFace& f2 = const_cast<VtkFace&>(f);
		f2.index = index++;
	}
	return ret;
}

struct ELineNotFound: public std::runtime_error{
	ELineNotFound(std::string s): std::runtime_error(s + " line not found while reading input") {};
};
std::string get_line_by_start(std::string start, std::istream& is){
	std::string line;
	while (is){
		std::getline(is, line);
		if (line.substr(0, start.size()) == start){
			return line;
		}
	}
	throw ELineNotFound(start);
}

std::pair<Point, double> face_normal_area(const std::vector<Point>& coo){
	if (coo.size() == 1){
		// point
		return { {1.0, 0.0, 0.0}, 1.0 };
	} else if (coo.size() == 2){
		// segment
		Vector v = coo[1] - coo[0];
		Vector ret = {v.y, -v.x, 0};
		double ret_len = vector_len(ret);
		return { ret/ret_len, ret_len };
	} else {
		// polygon
		Point ret;
		for (size_t i=1; i<coo.size()-1; ++i){
			Vector v1 = coo[i] - coo[0];
			Vector v2 = coo[i+1] - coo[0];
			ret += vector_cross(v1, v2);
		}
		double ret_len = vector_len(ret);
		return { ret/ret_len, ret_len/2 };
	}
}

Point aver_point(const std::vector<Point>& coo){
	Point ret;
	for (auto& p: coo) ret += p;
	ret /= coo.size();
	return ret;
}

double cell_volume(const std::vector<Point>& coo, CellCode code){
	switch (code){
		case CellCode::SEGMENT:
			return std::abs(coo[1].x - coo[0].x);
		case CellCode::POLYGON:
			return face_normal_area(coo).second;
		case CellCode::TETRAHEDRON:
			return vector_dot_triple(
				coo[1] - coo[0],
				coo[2] - coo[0],
				coo[3] - coo[0])/6;
		default:
			_THROW_NOT_IMP_;
	}
}

}


std::shared_ptr<UnstructuredGrid> UnstructuredGrid::read_from_vtk(std::string fn){
	std::cout << "Reading grid from " << fn << std::endl;
	
	std::ifstream ifs(fn);
	if (!ifs){
		throw std::runtime_error(fn + " is not found");
	}
	std::string line, tmp;
	// header
	line = get_line_by_start("DATASET", ifs);
	if (line.substr(8) != "UNSTRUCTURED_GRID"){
		throw std::runtime_error("Only unstructured grid can be read");
	}
	// points
	line = get_line_by_start("POINTS", ifs);
	int n_points;
	std::istringstream(line) >> tmp >> n_points;

	std::vector<Point> points(n_points);
	for (int i=0; i<n_points; ++i){
		ifs >> points[i].x >> points[i].y >> points[i].z;
	}
	std::cout << n_points << " points" << std::endl;

	// cells
	int n_totals, n_cells;
	line = get_line_by_start("CELLS", ifs);
	std::istringstream(line) >> tmp >> n_cells >> n_totals;
	std::vector<int> totals(n_totals);
	std::vector<int> types(n_cells);
	std::vector<int> bnd(n_cells, -1);
	// -- totals
	for (int i=0; i<n_totals; ++i){
		ifs >> totals[i];
	}
	// -- types
	line = get_line_by_start("CELL_TYPES", ifs);
	for (int i=0; i<n_cells; ++i){
		ifs >> types[i];
	}
	// -- bnd if exists
	try {
		line = get_line_by_start("SCALARS CellEntityIds", ifs);
		line = get_line_by_start("LOOKUP_TABLE", ifs);
		for (int i=0; i<n_cells; ++i){
			ifs >> bnd[i];
		}
	} catch (ELineNotFound){
	}

	// assemble vtk cells
	std::vector<std::shared_ptr<AVtkCell>> vtk_cells;
	int* cursor = &totals[0];
	int geometry_dimension = 0;
	for (int i=0; i<n_cells; ++i){
		int len = cursor[0];
		cursor++;
		vtk_cells.push_back(AVtkCell::build(
			std::vector<int>(cursor, cursor + len),
			types[i]));
		geometry_dimension = std::max(geometry_dimension, vtk_cells.back()->dim());
		cursor += len;
	}

	// get internal and boundary vtk cells
	std::vector<std::shared_ptr<AVtkCell>> internal_cells;
	std::map<int, std::vector<std::shared_ptr<AVtkCell>>> boundary_cells;
	int n_boundary_cells = 0;
	for (int i=0; i<n_cells; ++i){
		int dim = vtk_cells[i]->dim();
		if (dim == geometry_dimension){
			internal_cells.push_back(vtk_cells[i]);
		} else if (dim == geometry_dimension-1){
			boundary_cells[bnd[i]].push_back(vtk_cells[i]);
			n_boundary_cells += 1;
		}
	}
	std::cout << internal_cells.size() << " " << geometry_dimension << "D cells" << std::endl;

	// build faces
	VtkFaceData face_data = assemble_faces_vtk(internal_cells);
	
	// build return grid
	std::shared_ptr<UnstructuredGrid> ret(new UnstructuredGrid(geometry_dimension));

	// fill points
	ret->_points = points;

	// fill cells
	for (const auto& c: internal_cells) {
		ret->_cells.push_back(c->points);
		ret->_vtk_cell_codes.push_back(c->code());

		std::vector<Point> cellp;
		for (int ip: c->points){
			cellp.push_back(points[ip]);
		}
		double v = ::cell_volume(cellp, (CellCode)c->code());
		ret->_cell_volumes.push_back(v);
	}

	std::vector<int> internal_face_indices, boundary_face_indices;
	// fill faces and face->cell connectivity
	ret->_tab_face_cell.resize(face_data.n_faces());
	for (const VtkFace& f: face_data.faces){
		ret->_faces.push_back(f.points_vec());

		Vector normal;
		double area;
		std::vector<Point> facep;
		for (int ip: ret->_faces.back()){
			facep.push_back(points[ip]);
		}
		std::tie(normal, area) = face_normal_area(facep);
		ret->_face_normals.push_back(normal);
		ret->_face_areas.push_back(area);

		ret->_tab_face_cell[f.index].left_cell = f.left_cell;
		ret->_tab_face_cell[f.index].right_cell = f.right_cell;

		if (f.left_cell >= 0 && f.right_cell >= 0){
			internal_face_indices.push_back(ret->_faces.size());
		} else {
			boundary_face_indices.push_back(ret->_faces.size());
		}
	}

	ret->_internal_faces = internal_face_indices;
	ret->_boundary_faces = boundary_face_indices;

	return ret;
}

UnstructuredGrid::UnstructuredGrid(int dim): AGrid(dim){}

int UnstructuredGrid::n_points() const{
	return (int)_points.size();
}

Point UnstructuredGrid::point(int ipoint) const{
	return _points[ipoint];
}

int UnstructuredGrid::n_faces() const{
	return (int)_faces.size();
}

int UnstructuredGrid::n_cells() const{
	return (int)_cells.size();
}

int UnstructuredGrid::n_boundary_faces() const{
	return (int)_boundary_faces.size();
}

double UnstructuredGrid::face_area(int iface) const{
	return _face_areas[iface];
}

double UnstructuredGrid::cell_volume(int icell) const{
	return _cell_volumes[icell];
}

Vector UnstructuredGrid::face_normal(int iface) const{
	return _face_normals[iface];
}

std::vector<int> UnstructuredGrid::tab_face_point(int iface) const{
	return _faces[iface];
}

std::array<int, 2> UnstructuredGrid::tab_face_cell(int iface) const{
	return {_tab_face_cell[iface].left_cell, _tab_face_cell[iface].right_cell};
}

void UnstructuredGrid::vtk_save_cells(std::string fpath) const{
	std::ofstream ofs(fpath);

	// header
	ofs << "# vtk DataFile Version 2.0" << std::endl;
	ofs << "fvm output" << std::endl;
	ofs << "ASCII" << std::endl;

	// grid
	ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	ofs << "POINTS " << n_points() << " double" << std::endl;
	for (int ipoint=0; ipoint<n_points(); ++ipoint){
		auto point = this->point(ipoint);
		ofs << point.x << " " << point.y << " " << point.z << std::endl;
		std::cout << point.x << " " << point.y << " " << point.z << std::endl;
	}

	int nc = n_cells();

	std::vector<std::vector<int>> cell_array = vtk_cell_array();
	std::vector<int> cell_types = vtk_cell_types();

	int nfull = std::accumulate(
		cell_array.begin(), cell_array.end(), 0,
		[](int ret, const std::vector<int>& v)->int{ return ret + int(v.size()) + 1; });
	ofs << "CELLS " << nc << " " << nfull << std::endl;
	for (const auto& ca: cell_array){
		ofs << ca.size() << " ";
		for (int icell: ca){
			ofs << icell << " ";
		}
		ofs << std::endl;
	}

	ofs << "CELL_TYPES " << nc << std::endl;
	for (int v: cell_types){
		ofs << v << std::endl;
	}
}

void UnstructuredGrid::vtk_save_faces(std::string fpath) const{
	std::ofstream ofs(fpath);

	// header
	ofs << "# vtk DataFile Version 2.0" << std::endl;
	ofs << "fvm output" << std::endl;
	ofs << "ASCII" << std::endl;

	// grid
	ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	ofs << "POINTS " << n_points() << " double" << std::endl;
	for (int ipoint=0; ipoint<n_points(); ++ipoint){
		auto point = this->point(ipoint);
		ofs << point.x << " " << point.y << " " << point.z << std::endl;
		std::cout << point.x << " " << point.y << " " << point.z << std::endl;
	}

	int nc = n_faces();

	std::vector<std::vector<int>> cell_array;
	for (int iface=0; iface<n_faces(); ++iface){
		std::vector<int> ipoints = tab_face_point(iface);
		cell_array.push_back(ipoints);
	}

	int nfull = std::accumulate(
		cell_array.begin(), cell_array.end(), 0,
		[](int ret, const std::vector<int>& v)->int{ return ret + int(v.size()) + 1; });
	ofs << "CELLS " << nc << " " << nfull << std::endl;
	for (const auto& ca: cell_array){
		ofs << ca.size() << " ";
		for (int icell: ca){
			ofs << icell << " ";
		}
		ofs << std::endl;
	}

	int ct = -1;
	if (dim == 1){
		ct = 1;
	} else if (dim == 2){
		ct = 3;
	} else if (dim == 3){
		ct = 7;
	}
	ofs << "CELL_TYPES " << nc << std::endl;
	for (int i=0; i<nc; ++i){
		ofs << ct << std::endl;
	}
}

bool UnstructuredGrid::point_in_cell(Point p, int icell) const{
	_THROW_NOT_IMP_;
}

auto UnstructuredGrid::cell_centers_rtree() const -> rtree_t*{
	if (_cell_centers_rtree == nullptr){
		_cell_centers_rtree.reset(new rtree_t());
		for (int icell=0; icell<n_cells(); ++icell){
			Point cc = cell_center(icell);
			rtree_ret_t t {boost_point_t{cc.x, cc.y, cc.z}, icell};
			_cell_centers_rtree->insert(t);
		}
	}
	return _cell_centers_rtree.get();
}

int UnstructuredGrid::find_cell_index(Point p) const{
	constexpr int NRET = 5;

	auto* rtree = cell_centers_rtree();

	std::vector<rtree_ret_t> result_n;
	result_n.reserve(NRET);

	boost_point_t boost_p {p.x, p.y, p.z};

	rtree->query(bg::index::nearest(boost_p, NRET), std::back_inserter(result_n));

	for (auto& pi: result_n){
		if (point_in_cell(p, pi.second)){
			return pi.second;
		}
	}

	return -1;
}
