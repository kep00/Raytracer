#pragma once
#include "ray.h"
#include "rectangle.h"
#include "Sphere.h"
#include "../external/json.hpp"
#include "../external/simpleppm.h"
#include<string>
#include"arealight.h"
#include"pointlight.h"

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <random>

using json = nlohmann::json;
using namespace std;
using namespace Eigen;

class Output {
public:
	string filename;
	int size[2];
	Vector3f lookat;
	Vector3f up;
	Vector3f centre;
	Vector3f ai;
	Vector3f bkc;
	float fov, pt;

	int raypixel[2];
	bool global, antialiasing;
	int maxbounces;

};

class RayTracer {
public:

	vector<Sphere> s;
	vector<rectangle> r;
	vector <arealight> al;
	vector <pointlight> pl;
	//vector<light> s;
	vector <Output> op;
	
	


	void setValue(const json j);
	RayTracer(const json j);
	void run();
	
	Vector3f colorshpere(const obj& o, Vector3f& normal, Vector3f& inter, Ray& ray);
	Vector3f clamp(Vector3f v, float min, float max);
	bool isPointInShadow(Ray& ray, float dis);
	float random_double();
	bool isintersect(Ray& ray, obj& o, Vector3f& intersection_pt, Vector3f& normal);
	Vector3f rayCast(Ray& ray);
	Vector3f sample_hemisphere(const Vector3f& normal);
	Vector3f pathtrace(Ray& ray, int depth);
	
	
};