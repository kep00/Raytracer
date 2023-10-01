
#include <iostream>
#include "sphere.h"

using namespace std;
using namespace Eigen;

Sphere::Sphere(Eigen::Vector3f centre, double radius) {
    cout << "Created sphere with radius " << radius << " at (" <<centre<< ")." << endl;
    this->centre = centre;
    this->radius = radius;
}