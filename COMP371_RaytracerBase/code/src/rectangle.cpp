
#include "rectangle.h"



rectangle::rectangle(float p1a, float p1b, float p1c, float p2a, float p2b, float p2c, float p3a, float p3b, float p3c, float p4a, float p4b, float p4c) {
    this->p1[0] = p1a;
    this->p1[1] = p1b;
    this->p1[2] = p1c;
    this->p2[0] = p2a;
    this->p2[1] = p2b;
    this->p2[2] = p2c;
    this->p3[0] = p3a;
    this->p3[1] = p3b;
    this->p3[2] = p3c;
    this->p4[0] = p4a;
    this->p4[1] = p4b;
    this->p4[2] = p4c;
}

Vector3f rectangle::norm() {
	Eigen::Vector3f p21 = p2 - p1;
	Eigen::Vector3f p31 = p3 - p1;
	Eigen::Vector3f norm = p21.cross(p31);
	return norm;
}