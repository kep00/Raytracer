#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "obj.h"

using namespace std;
using namespace Eigen;

class Sphere: public obj
{
public:

	Eigen::Vector3f centre;
	float radius;

	Sphere(Eigen::Vector3f centre, double radius);
	

	
};