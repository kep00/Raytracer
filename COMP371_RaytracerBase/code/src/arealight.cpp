#include "arealight.h"




void arealight::setup(Vector3f p1, Vector3f p2, Vector3f p3, Vector3f p4, bool usecenter, int n) {
	this->p1 = p1;
	this->p2 = p2;
	this->p3 = p3;
	this->p4 = p4;
	this->usecenter = usecenter;
	this->n = n;

}

Vector3f arealight::centre() {
	float centerx = (p1[0] + p2[0] + p3[0] + p4[0]) / 4;
	float centery = (p1[1] + p2[1] + p3[1] + p4[1]) / 4;
	float centerz = (p1[2] + p2[2] + p3[2] + p4[2]) / 4;
	Vector3f centre(centerx, centery, centerz);
	return centre;
}