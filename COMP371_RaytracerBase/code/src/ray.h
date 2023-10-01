#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Ray {

    

public:
    Eigen::Vector3f orig;
    Eigen::Vector3f dir;
    Ray() {}
    Ray(const Vector3f origin, const Vector3f direction)
        : orig(origin), dir(direction)
    {}

    Vector3f origin() const { return orig; }
    Vector3f direction() const { return dir; }
    Vector3f at(float t) const {
        return orig + t * dir;
    }

   
};
