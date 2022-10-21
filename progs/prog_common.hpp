#ifndef CHECKS_HPP
#define CHECKS_HPP

#include "common.hpp"

// ============== default input directory
inline std::string from_input_path(std::string fname) {
	// CFDLIB_INPUT_DIR_PATH macro should be defined by the compiler
	return CFDLIB_INPUT_DIR_PATH + fname;
}

// ============== default output directory
inline std::string from_output_path(std::string fname) {
	// CFDLIB_OUTPUT_DIR_PATH macro should be defined by the compiler
	return CFDLIB_OUTPUT_DIR_PATH + fname;
}


// ============== checks for tests
bool equal_int_vec(const std::vector<int>& v1, const std::vector<int>& v2){
	if (v1.size() != v2.size()) return false;
	for (size_t i=0; i<v1.size(); ++i){
		if (v1[i] != v2[i]) return false;
	}
	return true;
}

bool equal_int_vec_anyorder(const std::vector<int>& v1, const std::vector<int>& v2){
	std::vector<int> v1s(v1);
	std::vector<int> v2s(v2);
	std::sort(v1s.begin(), v1s.end());
	std::sort(v2s.begin(), v2s.end());
	return equal_int_vec(v1s, v2s);
}

#define CHECK(cond)\
	{\
		if (!(cond)){ \
			printf("======= FAILED CHECK:\n"); \
			printf("function: %s\nat:       %s:%i\n", __PRETTY_FUNCTION__, __FILE__, __LINE__); \
			printf("failed:   %s\n", #cond);\
			throw std::runtime_error("check failed");\
		};\
	}\

#define CHECK_FLOAT(x, y, eps) CHECK(std::abs(x-y) < eps)

#define CHECK_FLOAT3(x, y) CHECK_FLOAT(x, y, 1e-3);

#endif
