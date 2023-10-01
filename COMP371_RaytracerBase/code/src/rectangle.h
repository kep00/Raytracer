#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "obj.h"

using namespace std;
using namespace Eigen;


class rectangle : public obj {
public:
	Vector3f p1,p2,p3,p4;

	rectangle(float p1a, float p1b, float p1c, float p2a, float p2b, float p2c, float p3a, float p3b, float p3c, float p4a, float p4b, float p4c);
	Vector3f norm();

};
