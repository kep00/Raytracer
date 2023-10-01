#include "pointlight.h"

using namespace std;
using Vector3f = Eigen::Vector3f;

pointlight::pointlight(Eigen::Vector3f& centre,bool& use) {
	this->centre = centre;
	this->use = use;

}