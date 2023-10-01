#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "light.h"

using namespace std;
using Vector3f = Eigen::Vector3f;

class pointlight : public Light {
public:

	Eigen::Vector3f centre;
	bool use;

	pointlight(Eigen::Vector3f& centre,bool& use);

};