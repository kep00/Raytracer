#pragma once

#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

class obj
{
public:
	Eigen::Vector3f ac, dc, sc;
	float ka, kd, ks;
	float pc;
	bool visible;
	string type;

	void setcolor(Eigen::Vector3f  ac, Eigen::Vector3f  dc, Eigen::Vector3f  sc);
	void setreflectioncoefficient(float & ka, float & kd, float & ks);
	void setPhongcoefficient(float & pc);
	



};
