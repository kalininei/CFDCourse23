#ifndef AELEMENT_HPP
#define AELEMENT_HPP

#include "common.hpp"
#include "geom.hpp"

class AElement{
public:
	~AElement() = default;
	AElement(const std::vector<int>& nbases, const std::vector<Point>& coo);

	virtual std::vector<double> mass() const;
	virtual std::vector<double> stiff() const;
	virtual std::vector<double> load() const;

	// total number of local bases
	int n_bases() const;

	// global index of local basis
	int i_basis(int ilocal) const;

	// total number of points in the cell
	int n_points() const;

	// shifted point coordinate. First is always 0
	Point point(int ilocal) const;
private:
	const std::vector<int> _bases;
	const std::vector<Point> _coo;
	const Point _point0;
};

#endif
