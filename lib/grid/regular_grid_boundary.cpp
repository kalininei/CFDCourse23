#include "grid/regular_grid_boundary.hpp"
#include "grid/regular_grid.hpp"

namespace{

std::vector<int> faces_by_direction(DirectionCode dircode, const ARegularGrid* grid){
	int nx = grid->nx();
	int ny = grid->ny();
	int nz = grid->nz();
	std::vector<int> ifaces;
	if (grid->dim == 1){
		switch (dircode){
			case DirectionCode::X_PLUS:
				ifaces.push_back(nx-1);
				break;
			case DirectionCode::X_MINUS:
				ifaces.push_back(0);
				break;
			default:
				_THROW_UNREACHABLE_;
		}
	} else if (grid->dim == 2){
		switch (dircode){
			case DirectionCode::X_PLUS:
				for (int iy=0; iy<ny-1; ++iy){
					ifaces.push_back(grid->face_ijk_to_glob(nx, iy, 0, 0));
				}
				break;
			case DirectionCode::X_MINUS:
				for (int iy=0; iy<ny-1; ++iy){
					ifaces.push_back(grid->face_ijk_to_glob(0, iy, 0, 0));
				}
				break;
			case DirectionCode::Y_PLUS:
				for (int ix=0; ix<nx-1; ++ix){
					ifaces.push_back(grid->face_ijk_to_glob(ix, ny, 0, 1));
				}
				break;
			case DirectionCode::Y_MINUS:
				for (int ix=0; ix<nx-1; ++ix){
					ifaces.push_back(grid->face_ijk_to_glob(ix, 0, 0, 1));
				}
				break;
			default:
				_THROW_UNREACHABLE_;
		}
	} else if (grid->dim == 3){
		switch (dircode){
			case DirectionCode::X_PLUS:
				for (int iy=0; iy<ny-1; ++iy)
				for (int iz=0; iz<nz-1; ++iz){
					ifaces.push_back(grid->face_ijk_to_glob(nx, iy, iz, 0));
				}
				break;
			case DirectionCode::X_MINUS:
				for (int iy=0; iy<ny-1; ++iy)
				for (int iz=0; iz<nz-1; ++iz){
					ifaces.push_back(grid->face_ijk_to_glob(0, iy, iz, 0));
				}
				break;
			case DirectionCode::Y_PLUS:
				for (int ix=0; ix<nx-1; ++ix)
				for (int iz=0; iz<nz-1; ++iz){
					ifaces.push_back(grid->face_ijk_to_glob(ix, ny, iz, 1));
				}
				break;
			case DirectionCode::Y_MINUS:
				for (int ix=0; ix<nx-1; ++ix)
				for (int iz=0; iz<nz-1; ++iz){
					ifaces.push_back(grid->face_ijk_to_glob(ix, 0, iz, 1));
				}
				break;
			case DirectionCode::Z_PLUS:
				for (int ix=0; ix<nx-1; ++ix)
				for (int iy=0; iy<ny-1; ++iy){
					ifaces.push_back(grid->face_ijk_to_glob(ix, iy, nz, 2));
				}
				break;
			case DirectionCode::Z_MINUS:
				for (int ix=0; ix<nx-1; ++ix)
				for (int iy=0; iy<ny-1; ++iy){
					ifaces.push_back(grid->face_ijk_to_glob(ix, iy, 0, 2));
				}
				break;
		}
	}

	return ifaces;
}
};

RegularGridBoundary::RegularGridBoundary(DirectionCode direction, const ARegularGrid* grid):
		AGridBoundary(faces_by_direction(direction, grid), grid), direction(direction){
}

RegularGridBoundary::RegularGridBoundary(DirectionCode direction, std::function<bool(Point)> filter, const ARegularGrid* grid):
		AGridBoundary(faces_by_direction(direction, grid), filter, grid), direction(direction){
}
