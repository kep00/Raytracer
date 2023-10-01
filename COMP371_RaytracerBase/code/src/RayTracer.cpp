#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <algorithm> 

#include "obj.h"
#include "Sphere.h"
#include "rectangle.h"
#include"RayTracer.h"
#include"light.h"
#include"pointlight.h"
#include"arealight.h"
#include "../external/json.hpp"
#include "../external/simpleppm.h"

#include <iomanip>
#define _USE_MATH_DEFINES // 使用math.h中的M_PI宏定义需要
#include <math.h>




using json = nlohmann::json;
using Vector3f = Eigen::Vector3f;
using namespace std;
std::random_device rd;
std::mt19937 gen(rd());





RayTracer::RayTracer(const json scene) {
    setValue(scene);

}


void RayTracer::setValue(const json j) {



    auto geometry = j["geometry"];
    auto light = j["light"];
    auto output = j["output"];

    // Access the values in the geometry array.
    for (const auto& shape : geometry) {
        std::string type = shape["type"];

        if (type == "rectangle") {


            auto ac = shape["ac"];
            auto dc = shape["dc"];
            auto sc = shape["sc"];

            float ka = shape["ka"];
            float kd = shape["kd"];
            float ks = shape["ks"];

            float pc = shape["pc"];


            Eigen::Vector3f sac(ac[0], ac[1], ac[2]);
            Eigen::Vector3f sdc(dc[0], dc[1], dc[2]);
            Eigen::Vector3f ssc(sc[0], sc[1], sc[2]);

            rectangle rectangle(shape["p1"][0], shape["p1"][1], shape["p1"][2], shape["p2"][0], shape["p2"][1], shape["p2"][2], shape["p3"][0], shape["p3"][1], shape["p3"][2], shape["p4"][0], shape["p4"][1], shape["p4"][2]);


            rectangle.setcolor(sac, sdc, ssc);
            rectangle.setreflectioncoefficient(ka, kd, ks);
            rectangle.setPhongcoefficient(pc);


            r.push_back(rectangle);

        }
        else if (type == "sphere") {
            auto centre = shape["centre"];
            auto radius = shape["radius"];

            auto ac = shape["ac"];
            auto dc = shape["dc"];
            auto sc = shape["sc"];

            float ka = shape["ka"];
            float kd = shape["kd"];
            float ks = shape["ks"];

            float pc = shape["pc"];
            Eigen::Vector3f scentre(centre[0], centre[1], centre[2]);
            Eigen::Vector3f sac(ac[0], ac[1], ac[2]);
            Eigen::Vector3f sdc(dc[0], dc[1], dc[2]);
            Eigen::Vector3f ssc(sc[0], sc[1], sc[2]);

            Sphere sphere(scentre, radius);
            sphere.setcolor(sac, sdc, ssc);
            sphere.setreflectioncoefficient(ka, kd, ks);
            sphere.setPhongcoefficient(pc);
            cout << "asads" << sphere.ac;
            s.push_back(sphere);

        }
    }

    // Access the values in the light array.
    for (const auto& l : light) {
        std::string type = l["type"];

        if (type == "point") {
            auto centre = l["centre"];
            auto id = l["id"];
            auto is = l["is"];
            bool use = true;
            if (l.contains("use")) {
                use = l["use"];
            }
            



            Eigen::Vector3f plcentre(centre[0], centre[1], centre[2]);
            Eigen::Vector3f plid(id[0], id[1], id[2]);
            Eigen::Vector3f plis(is[0], is[1], is[2]);
            pointlight pointlight(plcentre,use);
            pointlight.setid(plid);
            pointlight.setis(plis);

            pl.push_back(pointlight);


            // Do something with the point light data...
        }
        if (type == "area") {
            auto p1 = l["p1"];
            auto p2 = l["p2"];
            auto p3 = l["p3"];
            auto p4 = l["p4"];
            auto id = l["id"];
            auto is = l["is"];
            auto n = 0;
            auto usecenter = false;
            if (l.contains("usecenter")) {
                usecenter = l["usecenter"];
            }
            
            if (l.contains("n")) {
                n = l["n"];
            }
            

            Eigen::Vector3f p_1(p1[0], p1[1], p1[2]);
            Eigen::Vector3f p_2(p2[0], p2[1], p2[2]);
            Eigen::Vector3f p_3(p3[0], p3[1], p3[2]);
            Eigen::Vector3f p_4(p4[0], p4[1], p4[2]);
            Eigen::Vector3f alid(id[0], id[1], id[2]);
            Eigen::Vector3f alis(is[0], is[1], is[2]);
            arealight arealight;

            arealight.setup(p_1, p_2, p_3, p_4, usecenter, n);
            arealight.setid(alid);
            arealight.setis(alis);

            al.push_back(arealight);


        }
        std::cout << endl << pl.size() << "sdjiasjdioajjdwoij";

    }

    // Access the values in the output array.
    for (const auto& o : output) {
        std::string filename = o["filename"];
        auto size = o["size"];
        auto lookat = o["lookat"];
        auto up = o["up"];
        float fov = o["fov"];
        auto centre = o["centre"];
        auto ai = o["ai"];
        auto bkc = o["bkc"];

        Output outpu;
        outpu.filename = filename;

        outpu.size[0] = size[0];
        outpu.size[1] = size[1];

        outpu.lookat = Vector3f(lookat[0], lookat[1], lookat[2]);

        outpu.up = Vector3f(up[0], up[1], up[2]);
        outpu.fov = fov;
        outpu.centre = Vector3f(centre[0], centre[1], centre[2]);
        outpu.ai = Vector3f(ai[0], ai[1], ai[2]);
        outpu.bkc = Vector3f(bkc[0], bkc[1], bkc[2]);

        if (o.contains("raysperpixel")) {
            outpu.raypixel[0] = o["raysperpixel"][0];
            outpu.raypixel[1] = o["raysperpixel"][1];
        }
        if (o.contains("globalillum")) {
            outpu.global = o["globalillum"];
        }
        if (o.contains("maxbounces")) {
            outpu.maxbounces = o["maxbounces"];
        }
        if (o.contains("probterminate")) {
            outpu.pt = o["probterminate"];
        }
        if (o.contains("antialiasing")) {
            outpu.antialiasing = o["antialiasing"];
        }
        if (!o.contains("raysperpixel")) {
            outpu.raypixel[0] = 0;
            outpu.raypixel[1] = 0;
        }
        if (!o.contains("globalillum")) {
            outpu.global = false;
        }
        if (!o.contains("maxbounces")) {
            outpu.maxbounces = 0;
        }
        if (!o.contains("probterminate")) {
            outpu.pt = 0;
        }
        if (!o.contains("antialiasing")) {
            outpu.antialiasing = false;
        }
        op.push_back(outpu);


    }

}


void RayTracer::run() { 
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    std::mt19937 generator;
    cout << "Before in";
    cout << "Running ray tracer" << endl;
    for(auto& opt :op){
    int n = 1;
       if(al.size()>0){
    n = al[0].n;} 
    int a = opt.raypixel[0];
    int b = opt.raypixel[1];

    float width = opt.size[0];
    float height = opt.size[1];
    vector<double> buffer(3 * width * height);
   
    for (int y = height-1; y >=0 ; y--) {
        for (int x = 0; x < width; x++) {
            Vector3f finalcolor(0, 0, 0);
           
            if (opt.global ||opt.antialiasing) {
                for (int i = 0; i < n; ++i) {
                    for (int k = 0; k < a; ++k) {
                        for (int j = 0; j < b; ++j) {
                            float u = distribution(generator);
                            float v = distribution(generator);
                            Vector3f right = opt.lookat.cross(opt.up).normalized();
                            float delta = 2.0 * tan((opt.fov * M_PI / 180.0) / 2.0) / height;
                            Vector3f direction = opt.lookat + tan((opt.fov * M_PI / 180.0) / 2.0) * opt.up
                                - width / 2.0 * delta * right
                                + ((x + u) * delta) * right
                                - ((y + v) * delta) * opt.up;
                            direction = direction.normalized();
                            Ray ray(opt.centre, direction);
                            Vector3f color = rayCast(ray);
                            finalcolor += color;
                        }
                    }
                }
                finalcolor /= (n * a * b);
            }
            else {
                
                Vector3f right = opt.lookat.cross(opt.up).normalized();
                float delta = 2.0 * tan((opt.fov * M_PI / 180.0) / 2.0) / height;

                Vector3f direction = opt.lookat + tan((opt.fov * M_PI / 180.0) / 2.0) * opt.up - width / 2.0 * delta * right + (x * delta + delta / 2.0) * right - (y * delta + delta / 2.0) * opt.up;

                Ray ray(opt.centre, direction);
                finalcolor = rayCast(ray);


            }
            buffer[3 * y * width + 3 * x + 0] = finalcolor[0];
            buffer[3 * y * width + 3 * x + 1] = finalcolor[1];
            buffer[3 * y * width + 3 * x + 2] = finalcolor[2];
        }
    }

    save_ppm(opt.filename, buffer, width, height);
    }
}



bool hit_sphere(const Vector3f center, double radius, const Ray& r, float& t) {
    Vector3f oc = r.origin() - center;
    auto a = r.direction().dot(r.direction());
    auto b = 2.0 * (oc.dot(r.direction()));
    auto c = oc.dot(oc) - radius * radius;
    auto discriminant = b * b - 4 * a * c;
    float t1 = (-b - sqrt(discriminant)) / (2.0 * a);
    float t2 = (-b + sqrt(discriminant)) / (2.0 * a);

    if (t1 < 0) {
        t = t2;
    }
    else if (t2 < 0) {
        t = t1;
    }
    else {
        t = min(t1, t2);
    }

    return (discriminant > 0);
}


bool hit_rectangle(const Vector3f p1, const Vector3f p2, const Vector3f p3, const Vector3f p4, const Ray& r, float& t) {
    Vector3f normal = (p2 - p1).cross(p4 - p1).normalized();
    float d = -(normal.dot(p1));
    float denominator = normal.dot(r.direction());
    if (denominator == 0.0f) {
        return false;
    }
    float numerator = -(normal.dot(r.origin()) + d);
    float intersectionDistance = numerator / denominator;
    if (intersectionDistance < 0.0f) {
        return false;
    }
    Vector3f intersectionPoint = r.origin() + intersectionDistance * r.direction();
    // Check if the intersection point is inside the rectangle
    Vector3f v1 = p2 - p1;
    Vector3f v2 = p3 - p2;
    Vector3f v3 = p4 - p3;
    Vector3f v4 = p1 - p4;

    Vector3f c1 = intersectionPoint - p1;
    Vector3f c2 = intersectionPoint - p2;
    Vector3f c3 = intersectionPoint - p3;
    Vector3f c4 = intersectionPoint - p4;

    if (normal.dot(v1.cross(c1)) < 0.0f ||
        normal.dot(v2.cross(c2)) < 0.0f ||
        normal.dot(v3.cross(c3)) < 0.0f ||
        normal.dot(v4.cross(c4)) < 0.0f) {
        return false;
    }

    // Calculate the distance from the ray origin to the intersection point
    t = intersectionDistance;
    return true;

}

Vector3f RayTracer::colorshpere(const obj& o, Vector3f& normal, Vector3f& inter, Ray& ray) {
    Vector3f clamped(0, 0, 0);
    for (auto& opt : op) {
        Vector3f color, color1(0, 0, 0), color2(0, 0, 0);
        Vector3f  a;
        if (!opt.global) {
            a = o.ka * o.ac;
        }
        else
        {
            a = color1;
        }


        for (pointlight pointlight : pl) {
            if (pointlight.use) {
                Vector3f light = (pointlight.centre - inter).normalized();
                Ray r(inter + light * 1e-3, light);
                float dis = (pointlight.centre - inter).norm();
                if (!isPointInShadow(r, dis)) {
                    color1 += o.kd * (o.dc.cwiseProduct(pointlight.id * (normal.dot(light))));
                    Vector3f view = (ray.orig - inter).normalized();
                    Vector3f H = (light + view).normalized();

                    color1 += o.ks * (o.sc.cwiseProduct(pointlight.is * pow(normal.dot(H), o.pc)));

                }
            }
        }


        for (arealight pointlight : al) {

            if (opt.global || pointlight.usecenter) {
                Vector3f light = (pointlight.centre() - inter).normalized();
                Ray r(inter + light * 1e-3, light);
                float dis = (pointlight.centre() - inter).norm();
                if (!isPointInShadow(r, dis)) {
                    color2 += o.kd * (o.dc.cwiseProduct(pointlight.id * (normal.dot(light))));

                    if (!opt.global) {
                        Vector3f view = (ray.orig - inter).normalized();
                        Vector3f H = (light + view).normalized();

                        color2 += o.ks * (o.sc.cwiseProduct(pointlight.is * pow(normal.dot(H), o.pc)));
                    }
                }
            }
            else
            {
                int a = opt.raypixel[0];
                int b = opt.raypixel[1];

                for (int i = 0; i < pointlight.n; ++i) {
                    for (int k = 0; k < a; ++k) {
                        for (int j = 0; j < b; ++j) {
                            std::uniform_real_distribution<float> dist(0, 1);
                            float u = dist(gen);
                            float v = dist(gen);
                            Vector3f edge1 = pointlight.p2 - pointlight.p1;
                            Vector3f edge2 = pointlight.p4 - pointlight.p1;
                            Vector3f centre = pointlight.p1 + u * edge1 + v * edge2;

                            Vector3f light = (centre - inter).normalized();
                            Ray r(inter + light * 1e-3, light);
                            float dis = (centre - inter).norm();
                            if (!isPointInShadow(r, dis)) {

                                float lambert = std::max(0.0f, normal.dot(light));
                                color2 += (o.kd * (o.dc.cwiseProduct(pointlight.id * (normal.dot(light)))));
                                Vector3f view = (ray.orig - inter).normalized();
                                Vector3f H = (light + view).normalized();

                                color2 += o.ks * (o.sc.cwiseProduct(pointlight.is * pow(normal.dot(H), o.pc)));
                            }
                        }
                    }
                }
                color2 = color2 / static_cast<float>(pointlight.n * a * b);
            }



        }

        color = a + color1 + color2;
         clamped = clamp(color, 0.0, 1.0); 
        return clamped;
    }

    return clamped;
}


Vector3f RayTracer::clamp(Vector3f v, float min, float max) {
    for (int i = 0; i < 3; i++) {
        v[i] = std::min(std::max(v[i], min), max);
    }
    return v;
}

Vector3f RayTracer::rayCast(Ray& ray) {
    Vector3f color(0, 0, 0);
    for (auto& opt : op) {
        Vector3f intersection_point(0, 0, 0);
        obj closest_shape;
        bool intersectionFound = false;
        Vector3f normal(0, 0, 0);
        
        if (opt.global) {
            color = pathtrace(ray, 0);
            color = clamp(color, 0.0, 1.0);
        }
        else
        {


            intersectionFound = isintersect(ray, closest_shape, intersection_point, normal);
            if (intersectionFound) {
                color = colorshpere(closest_shape, normal, intersection_point, ray);
            }
            else {
                color = opt.bkc;
            }

        }
        return color;
    }
    return color;
}
bool RayTracer::isintersect(Ray& ray, obj& o, Vector3f& intersection_pt, Vector3f& normal) {
    float t_cur = numeric_limits<float>::max();
    bool intersect = false;
    for (const auto& sphere : s) {
        float t;
        if (hit_sphere(sphere.centre, sphere.radius, ray, t) && t < t_cur) {
            t_cur = t;
            intersect = true;
            o = sphere;
            normal = (ray.at(t_cur) - sphere.centre).normalized();
        }
    }

    for (auto& rectangle : r) {
        float t = 0;
        if (hit_rectangle(rectangle.p1, rectangle.p2, rectangle.p3, rectangle.p4, ray, t) && t < t_cur) {
            t_cur = t;
            Vector3f inter = ray.at(t_cur);
            intersect = true;
            o = rectangle;
            normal = rectangle.norm().normalized();

        }

    }

    if (intersect) {
        intersection_pt = ray.at(t_cur);

    }
    return intersect;
}


bool RayTracer::isPointInShadow(Ray& ray, float dis) {
    float t = 0.0;
    for (const auto& sphere : s) {
        if (hit_sphere(sphere.centre, sphere.radius, ray, t))
        {
            if (t > 0.00001 && t < dis) return true;
        }

    }
    for (auto& rectangle : r) {
        if (hit_rectangle(rectangle.p1, rectangle.p2, rectangle.p3, rectangle.p4, ray, t)) {
            if (t > 0.00001 && t < dis) return true;
        }
    }
    return false;

}



float RayTracer::random_double() {
    std::uniform_real_distribution<float> dist(0, 1);
    return dist(gen);


}

Vector3f RayTracer::sample_hemisphere(const Vector3f& normal) {

    float u1 = random_double();
    float u2 = random_double();
    float r = sqrt(u1);
    float theta = 2.0 * M_PI * u2;

    Vector3f sample(r * cos(theta), r * sin(theta), sqrt(1.0 - u1));

    Vector3f up = std::abs(normal.z()) < 0.9 ? Vector3f::UnitZ() : Vector3f::UnitX();
    Vector3f tangent = up.cross(normal).normalized();
    Vector3f bitangent = normal.cross(tangent);

    return tangent * sample.x() + bitangent * sample.y() + normal * sample.z();
}

Vector3f RayTracer::pathtrace(Ray& ray, int depth) {
    Vector3f indirectIllumination(0, 0, 0);
    Vector3f directIllumination(0, 0, 0);
    for (auto& opt : op) {
        if (depth >= opt.maxbounces)
        {
            return Vector3f(0, 0, 0);
        }
        obj o;
        Vector3f inter(0, 0, 0), normal(0, 0, 0);
        if (!isintersect(ray, o, inter, normal)) {
            return opt.bkc;
        }
        Vector3f randomDirection = sample_hemisphere(normal);
        float cos_theta = std::max(normal.dot(randomDirection), 0.0f);
        
        directIllumination = colorshpere(o, normal, inter, ray);
        if (random_double() < opt.pt) {

            return directIllumination;
        }


        Vector3f brdf = o.dc * o.kd / M_PI;
        Ray newray(inter + randomDirection * 1e-3, randomDirection);

        indirectIllumination += ((brdf.cwiseProduct(pathtrace(newray, depth + 1)) * cos_theta * (2 * M_PI)));

        return (directIllumination + indirectIllumination);
    }
    return (directIllumination + indirectIllumination);
}