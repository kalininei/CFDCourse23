#ifndef _FE_QUADRATURE_HPP_
#define _FE_QUADRATURE_HPP_

#include <array>
#include <functional>
#include "geom.hpp"

class IQuadratureRule{
public:
	virtual ~IQuadratureRule() = default;
	virtual double integrate(const std::function<double(Point)>& func) const = 0;
	virtual std::vector<double> integrate(const std::function<std::vector<double>(Point)>& func) const = 0;
};

template<size_t N>
class AQuadratureRule: public IQuadratureRule{
public:
	virtual ~AQuadratureRule() = default;

	double integrate(const std::function<double(Point)>& func) const override{
		double ret = 0;
		for (size_t i=0; i<N; ++i){
			ret += this->weights[i] * func(coo[i]);
		}
		return ret;
	}
	std::vector<double> integrate(const std::function<std::vector<double>(Point)>& func) const override{
		std::vector<double> ret = func(coo[0]);
		for (double& k: ret) k *= weights[0];
		for (size_t i=1; i<N; ++i){
			std::vector<double> v = func(coo[i]);
			for (size_t j=0; j<ret.size(); ++j){
				ret[j] += this->weights[i] * v[j];
			}
		}
		return ret;
	}
protected:
	AQuadratureRule(std::array<double, N>&& w, std::array<Point, N>&& coo);
private:
	std::array<double, N> weights;
	std::array<Point, N> coo;
};


// ================= [-1, 1] segment
class SegmentQuad1: public AQuadratureRule<1>{
public:
	SegmentQuad1();
};
class SegmentQuad2: public AQuadratureRule<2>{
public:
	SegmentQuad2();
};
class SegmentQuad3: public AQuadratureRule<3>{
public:
	SegmentQuad3();
};
class SegmentQuad4: public AQuadratureRule<4>{
public:
	SegmentQuad4();
};

// ================= [-1, 1]^2 square
class SquareQuad1: public AQuadratureRule<1>{
public:
	SquareQuad1();
};
class SquareQuad2: public AQuadratureRule<4>{
public:
	SquareQuad2();
};
class SquareQuad3: public AQuadratureRule<9>{
public:
	SquareQuad3();
};
class SquareQuad4: public AQuadratureRule<16>{
public:
	SquareQuad4();
};

// ================= [-1, 1]^3 cube
class CubeQuad1: public AQuadratureRule<1>{
public:
	CubeQuad1();
};
class CubeQuad2: public AQuadratureRule<8>{
public:
	CubeQuad2();
};
class CubeQuad3: public AQuadratureRule<27>{
public:
	CubeQuad3();
};
class CubeQuad4: public AQuadratureRule<64>{
public:
	CubeQuad4();
};

// ================= [(0, 0)-(1, 0)-(0, 1)] triangle
class TriQuad1: public AQuadratureRule<1>{
public:
	TriQuad1();
};
class TriQuad2: public AQuadratureRule<3>{
public:
	TriQuad2();
};
class TriQuad3: public AQuadratureRule<4>{
public:
	TriQuad3();
};
class TriQuad4: public AQuadratureRule<6>{
public:
	TriQuad4();
};

// ================= [(0, 0, 0)-(1, 0, 0)-(0, 1, 0-(0, 0, 1)] tetrahedron


#endif
