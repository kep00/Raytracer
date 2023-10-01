#include "obj.h"

using namespace std;
using namespace Eigen;

void obj::setcolor(Eigen::Vector3f ac, Eigen::Vector3f dc, Eigen::Vector3f sc) {
	this->ac = ac;
	this->dc = dc;
	this->sc = sc;
}

void obj::setreflectioncoefficient(float &ka, float &kd, float &ks) {
	this->ka = ka;
	this->kd = kd;
	this->ks = ks;

}
void obj::setPhongcoefficient(float &pc) {

	this->pc = pc;
}