#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Light {
public:
	Vector3f id;
	Vector3f is;

	void setid(Vector3f id) { this->id = id; };
	void setis(Vector3f is) { this->is = is; };

};