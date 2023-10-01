#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>
#include "light.h"
#include <vector>
#include <random>

using namespace std;
using namespace Eigen;


class arealight : public Light {
public:
	Eigen::Vector3f p1, p2, p3, p4;
    bool usecenter;
    int n;
    
	void setup(Vector3f p1, Vector3f p2, Vector3f p3, Vector3f p4,bool usecenter, int n);
   
	Vector3f centre();
	Vector3f randomPoint() {
		float u = ((float)rand() / (RAND_MAX));
		float v = ((float)rand() / (RAND_MAX));
		Vector3f random_point = centre() + u * p1 + v * p2;
		return random_point;
	}
    
   
    Vector3f sample(int num_samples)  {
       
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dist_width(-(p1[0] - p2[0]) / 2, (p1[0] - p2[0]) / 2);
        std::uniform_real_distribution<float> dist_height(-(p1[1] - p2[1]) / 2, (p1[1] - p2[1]) / 2);

        
            float x = dist_width(gen);
            float y = dist_height(gen);
            Vector3f a(x, y, 0);
            Vector3f sample_pos = centre() + a;
            

        return sample_pos;
    }

};