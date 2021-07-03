#define _USE_MATH_DEFINES
#include <limits>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include "geometry.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

int envmap_width, envmap_height;
std::vector<Vec3f> envmap;
const float noise_amplitude = 0.4;
bool ph_ray_intersect(Vec3f orig, Vec3f dir, float& t0) {
    float radius = 100.0f;
    Vec3f center = Vec3f (0,0,0);
    Vec3f L = center - orig;
    float tca = L*dir;
    float d2 = L*L - tca*tca;
    if (d2 > radius*radius) return false;
    float thc = sqrtf(radius*radius - d2);
    t0       = tca - thc;
    float t1 = tca + thc;
    if (t0 < 0) t0 = t1;
    if (t0 < 0) return false;
    return true;
}
template <typename T> inline T lerp(const T &v0, const T &v1, float t) {
    return v0 + (v1-v0)*std::max(0.f, std::min(1.f, t));
}
Vec3f palette_fire(const float d) {
    const Vec3f     white(1.2, 1.2, 1.2);
    const Vec3f   yellow(0.9, 0.9, 0.2); // note that the color is "hot", i.e. has components >1
    const Vec3f   orange(0.7, 0.42, 0.0);
    const Vec3f      red(0.7, 0.0, 0.0);
    const Vec3f darkgray(0.15, 0.15, 0.15);
    const Vec3f     gray(0.3, 0.3, 0.3);
    float x = std::max(0.f, std::min(2.5f, d));
    if (x<.25f)
        return lerp(gray, darkgray, x*4.f)*1.2;
    else if (x<.5f)
        return lerp(darkgray, red, x*4.f-1.f)*1.2;
    else if (x<1.f)
        return lerp(red, orange, x*2.f-1.f)*1.2;
    else if (x<1.75f)
        return lerp(orange, yellow, x*1.33f-1.33f)*1.2;

    return lerp(yellow, white, x*4.f-7.f)*1.2;
}
float hash(const float n) {
    float x = sin(n)*43758.5453f;
    return x-floor(x);
}

float noise(const Vec3f &x) {
    Vec3f p(floor(x.x), floor(x.y), floor(x.z));
    Vec3f f(x.x-p.x, x.y-p.y, x.z-p.z);
    f = f*(f*(Vec3f(3.f, 3.f, 3.f)-f*2.f));
    float n = p*Vec3f(1.f, 57.f, 113.f);
    return lerp(lerp(
            lerp(hash(n +  0.f), hash(n +  1.f), f.x),
            lerp(hash(n + 57.f), hash(n + 58.f), f.x), f.y),
                lerp(
                        lerp(hash(n + 113.f), hash(n + 114.f), f.x),
                        lerp(hash(n + 170.f), hash(n + 171.f), f.x), f.y), f.z);
}

Vec3f rotate(const Vec3f &v) {
    return Vec3f(Vec3f(0.00,  0.80,  0.60)*v, Vec3f(-0.80,  0.36, -0.48)*v, Vec3f(-0.60, -0.48,  0.64)*v);
}

float fractal_brownian_motion(const Vec3f &x) {
    Vec3f p = rotate(x);
    float f = 0;
    f += 0.5000*noise(p); p = p*2.32;
    f += 0.2500*noise(p); p = p*3.03;
    f += 0.1250*noise(p); p = p*2.61;
    f += 0.0625*noise(p);
    return f/0.9375;
}
struct Material {
    Material(const float &r, const Vec4f &a, const Vec3f &color, const float &spec) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : refractive_index(1), albedo(1,0,0,0), diffuse_color(), specular_exponent() {}
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};
struct Light {
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};
Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

Vec3f refraction(const Vec3f &I, Vec3f N,float n1, float n2, bool& all_reflect){
    if (I*N<0.0f){
        N = -N;
    }
    float cosA = I*N;
    float sinA = sqrt(1-cosA*cosA);
    Vec3f q = (I-N*cosA);
    Vec3f tang = q.normalize();
    float sinB = sinA/ n2 * n1;
    all_reflect= false;
    if (sinB>0.999){
        all_reflect = true;
        sinB = 0.999;
    }
    if (sinB<-0.999){
        all_reflect = true;
        sinB = -0.999;
    }
    float cosB = sqrt(1-sinB*sinB);
    return  tang * sinB + N*cosB ;
}
float help_signed_distance(Vec3f p) {
    Vec3f v = p;
    v.z=v.z+3;//координаты центра (0,0.75,-3)
    v.y = v.y ;
    v.x = v.x;
    float displacement =0;// -fractal_brownian_motion(p*3.4)*noise_amplitude;;
    return v.norm() - (0.8f + displacement);
}
float boom_signed_distance(Vec3f p) {
    Vec3f v = p;
    v.z=v.z+3;//координаты центра (0,0.75,-3)
    v.y = v.y ;
    v.x = v.x;
    float displacement =fractal_brownian_motion(p*3.4)*noise_amplitude;;
    return v.norm() - (0.4f + displacement);
}
Vec3f boom_field_normal(Vec3f &pos) {
    const float eps = 0.01;
    float d = boom_signed_distance(pos);
    float nx = boom_signed_distance(pos + Vec3f(eps, 0, 0)) - d;
    float ny = boom_signed_distance(pos + Vec3f(0, eps, 0)) - d;
    float nz = boom_signed_distance(pos + Vec3f(0, 0, eps)) - d;
    return Vec3f(nx, ny, nz).normalize();
}
float torus_signed_distance(Vec3f p) {
    Vec3f v = p;
    v.z=v.z+3;//координаты центра (0,0,-3)
    v.x = v.x;
    v.y = v.y+1.5;

    Vec2f t(1,0.25);
    Vec3f b(v.z,v.x,0);
    Vec3f q = Vec3f(b.norm()-t.x,v.y,0);
    return q.norm()-t.y;

}
Vec3f torus_field_normal(Vec3f &pos) {
    const float eps = 0.01;
    Vec3f z1 = pos + Vec3f(eps, 0, 0);
    Vec3f z2 = pos - Vec3f(eps, 0, 0);
    Vec3f z3 = pos + Vec3f(0, eps, 0);
    Vec3f z4 = pos - Vec3f(0, eps, 0);
    Vec3f z5 = pos + Vec3f(0, 0, eps);
    Vec3f z6 = pos - Vec3f(0, 0, eps);
    float dx = torus_signed_distance(z1) - torus_signed_distance(z2);
    float dy = torus_signed_distance(z3) - torus_signed_distance(z4);
    float dz = torus_signed_distance(z5) - torus_signed_distance(z6);
    return Vec3f(dx, dy, dz).normalize();
}
float cube_signed_distance(Vec3f p) {
    Vec3f v = p;
    v.z=abs(v.z+5)/5;//координаты центра (0,0,-3)
    v.y =abs( v.y+2.75);
    v.x = abs(v.x)/10;
    return std::max(std::max(v.x,v.y),v.z)-1;
}
Vec3f cube_field_normal(Vec3f &pos) {
    const float eps = 0.1;
    Vec3f z1 = pos + Vec3f(eps, 0, 0);
    Vec3f z2 = pos - Vec3f(eps, 0, 0);
    Vec3f z3 = pos + Vec3f(0, eps, 0);
    Vec3f z4 = pos - Vec3f(0, eps, 0);
    Vec3f z5 = pos + Vec3f(0, 0, eps);
    Vec3f z6 = pos - Vec3f(0, 0, eps);
    float dx = cube_signed_distance(z1) - cube_signed_distance(z2);
    float dy = cube_signed_distance(z3) - cube_signed_distance(z4);
    float dz = cube_signed_distance(z5) - cube_signed_distance(z6);
    return (Vec3f(dx/ (2.0*eps), dy/ (2.0*eps), dz/ (2.0*eps)) ).normalize();
}
float stand_signed_distance(Vec3f p) {
    Vec3f v = p;
    v.z=v.z+6;//координаты центра (0,0,-3)
    v.y=(v.y+1.75)*4;
    v.x=(v.x+3);
    Vec3f q(v.x,v.z,0);
    float h = 1;
    float r = 1;
    Vec3f g = Vec3f (q.norm(),abs(v.y),0)-Vec3f(h,r,0);
    float t = Vec3f(std::max(g.x,0.0f),std::max(g.y,0.0f),0).norm();
    return t+std::min(std::max(g.x,g.y),0.0f) ;
}
Vec3f stand_field_normal(Vec3f &pos) {
    const float eps = 0.01;
    Vec3f z1 = pos + Vec3f(eps, 0, 0);
    Vec3f z2 = pos - Vec3f(eps, 0, 0);
    Vec3f z3 = pos + Vec3f(0, eps, 0);
    Vec3f z4 = pos - Vec3f(0, eps, 0);
    Vec3f z5 = pos + Vec3f(0, 0, eps);
    Vec3f z6 = pos - Vec3f(0, 0, eps);
    float dx = stand_signed_distance(z1) - stand_signed_distance(z2);
    float dy = stand_signed_distance(z3) - stand_signed_distance(z4);
    float dz = stand_signed_distance(z5) - stand_signed_distance(z6);
    return (Vec3f(dx/ (2.0*eps), dy/ (2.0*eps), dz/ (2.0*eps)) ).normalize();
}
float glass_signed_distance(Vec3f p) {
    Vec3f v = p;
    v.z=v.z+3;//координаты центра (0,0,-3)
    v.y = v.y ;
    v.x = v.x;
    return v.norm() - 1.55;
}
Vec3f glass_field_normal(Vec3f &pos) {
    const float eps = 0.01;
    float d = glass_signed_distance(pos);
    float nx = glass_signed_distance(pos + Vec3f(eps, 0, 0)) - d;
    float ny = glass_signed_distance(pos + Vec3f(0, eps, 0)) - d;
    float nz = glass_signed_distance(pos + Vec3f(0, 0, eps)) - d;
    return Vec3f(nx, ny, nz).normalize();
}
bool scene_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N, Material &material, int in_mat, int& out_mat, Vec3f& pp_color, int ignore=-1) {
    float min_dist = std::numeric_limits<float>::max();
    hit = orig; // TODO
    int n_obj = 5;
    std::vector<float> dists(n_obj,min_dist);
    int id_first_hit=-1;
    float d=100000.0f;
    for (size_t i=0; i<512; i++) {
        dists[0] = help_signed_distance(hit);
        dists[1] = torus_signed_distance(hit);
        dists[2] = cube_signed_distance(hit);
        dists[3] = stand_signed_distance(hit);
        dists[4] = glass_signed_distance(hit);
        if (ignore!=-1){
            dists[ignore]= 100000.0f;
        }
        if (in_mat!=-1){
            dists[in_mat]*=-1;
        }
        for (int q=0;q<n_obj;q++){
            if (dists[q]<d){
                id_first_hit = q;
                d = dists[q];
            }
        }
        if (id_first_hit==0){
            d = boom_signed_distance(hit);
        }

        out_mat = id_first_hit;
        if (out_mat == in_mat){
            out_mat = -1;
        }
        if (in_mat==0){
            out_mat =0;
            while(dists[0]>0){
                d = 0.05;
                Vec3f v = hit;
                hit = hit + dir*std::max(d,0.01f);
                dists[0] = -1*boom_signed_distance(hit);
                v.z= v.z+3;//координаты центра (0,0.75,-3)
                v.y = v.y ;
                v.x = v.x;
                float noise_level = (0.8f-v.norm())/noise_amplitude;
                Vec3f curr_color;
                curr_color = palette_fire((-.2 + noise_level)*2);
                pp_color  = pp_color*0.8+curr_color*0.2;
                out_mat++;
            }
            return true;
            if (dists[0]>0){
                d = 0.1;
                Vec3f curr_color;
                Vec3f v = hit;
                v.z= v.z+3;//координаты центра (0,0.75,-3)
                v.y = v.y ;
                v.x = v.x;
                float noise_level = (0.8f-v.norm())/noise_amplitude;
                curr_color = palette_fire((-.2 + noise_level)*2);
                pp_color  = curr_color;
                //pp_color.x = pp_color.x*pow(curr_color.x,0.1);
                //pp_color.y = pp_color.y*pow(curr_color.y,0.1);
                //pp_color.z = pp_color.z*pow(curr_color.z,0.1);
                //накопление цвета от огня
            }
        }
        hit = hit + dir*std::max(d,0.01f);
        if (d<0.001){
            switch (id_first_hit) {
                case 0: // Взрыв
                    material.albedo = Vec4f (0.9,0.1,0,0);
                    material.specular_exponent = 100;
                    material.diffuse_color = Vec3f (1,1,1);
                    material.refractive_index = 1.;
                    N = boom_field_normal(hit);
                    break;
                case 1: // Тор
                    material.albedo = Vec4f (0.9,0.6,0.3,0);
                    material.specular_exponent = 100;
                    material.diffuse_color = Vec3f (0.2,0.05,0.05);
                    material.refractive_index = 1.;
                    N = torus_field_normal(hit);
                    break;
                case 2: // Куб
                    material.albedo = Vec4f (0.9,0.1,0.1,0);
                    material.specular_exponent = 100;
                    material.diffuse_color = Vec3f (0.2,0.2,0.8);
                    material.refractive_index = 1.;

                    N = cube_field_normal(hit);
                    break;
                case 3: // Цилиндр
                    material.albedo = Vec4f (0.9,0.6,0.3,0);
                    material.specular_exponent = 100;
                    material.diffuse_color = Vec3f (0.85,0.3,0.3);
                    material.refractive_index = 1.;

                    N = stand_field_normal(hit);
                    break;
                case 4: // Стекло
                    material.albedo = Vec4f(0.0,  1.6, 0.1, 0.6);
                    material.specular_exponent = 150.;
                    material.diffuse_color = Vec3f(0.6, 0.7, 0.8);
                    material.refractive_index = 1.5;

                    N = glass_field_normal(hit);
                    break;
                default: // если никуда не попали
                    break;
            }
            return true;
        }
    }
    hit = orig;
    return false;

}
Vec3f ray_marching(const Vec3f &orig, const Vec3f &dir,int cur_mat, size_t depth=0){
    Vec3f point, N;
    Material material;
    Vec3f tmp_vec;
    int out_m = -1;

    if (depth>4 || !scene_intersect(orig, dir, point, N, material,cur_mat,out_m,tmp_vec)) {
        if (depth>4){
            return Vec3f(0.4, 0.4, 0.4); // background color

        }
        float dist = 0;
        ph_ray_intersect(orig, dir, dist);
        Vec3f p = orig+dir*dist;
        int a = (atan2(p.z, p.x)/(2*M_PI) + .5)*envmap_width;
        int b = acos(p.y/100)/M_PI*envmap_height;
        return envmap[a+b*envmap_width];
        //TODO
    }
    if (out_m==0){
        Vec3f v = point;
        Vec3f new_point;
        v.z= v.z+3;//координаты центра (0,0.75,-3)
        v.y = v.y ;
        v.x = v.x;
        float noise_level = (0.8f-v.norm())/noise_amplitude;
        int tmp_mat=4;
        Vec3f pp_color(1,1,1);
        scene_intersect(point+dir*0.01, dir, new_point, N, material,out_m,tmp_mat,pp_color);
        //return palette_fire((-.2 + noise_level)*2);
        //Vec3f res=ray_marching(orig, dir, out_m,depth + 1);
        Vec3f res = Vec3f(1.1, 1, 1)*1.2;
        float x = pow(0.8,tmp_mat);
        res = Vec3f (res.x*pp_color.x,res.y*pp_color.y,res.z*pp_color.z)*(1-x)+ray_marching(new_point, dir, 4,depth + 1)*x;
        return res;
    }
    float n1 = 1.5;
    float n2 = 1.0;
    Vec3f refract_dir;
    Vec3f reflect_dir;
    Vec3f refract_orig;
    Vec3f refract_color;
    if (cur_mat ==-1){
        bool all_ref_flag;
        refract_dir = refraction(dir, N, n2,n1,all_ref_flag).normalize();
        refract_orig = refract_dir*N < 0 ? point - N*5e-2 : point + N*5e-2;
        refract_color = ray_marching(refract_orig, refract_dir, out_m,depth + 1);
    }else{
        bool all_ref_flag;
        refract_dir = refraction(dir, N, n1,n2,all_ref_flag).normalize();
        if (all_ref_flag){
            refract_dir = reflect(dir, N).normalize();
            refract_orig = refract_dir*N < 0 ? point - N*5e-2 : point + N*5e-2;
            refract_color = ray_marching(refract_orig, refract_dir, cur_mat,depth + 1);
        }else{
            refract_orig = refract_dir*N < 0 ? point - N*5e-2 : point + N*5e-2;
            refract_color = ray_marching(refract_orig, refract_dir, out_m,depth + 1);
        }
    }
    reflect_dir = reflect(dir, N).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*5e-2 : point + N*5e-2; // offset the original point to avoid occlusion by the object itself
    Vec3f reflect_color = ray_marching(reflect_orig, reflect_dir, cur_mat,depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    std::vector<Light> lights;
    lights.push_back(Light(Vec3f (-21, 7, 12),0.3));
    lights.push_back(Light(Vec3f (12, 5, 2),0.3));

    for (int i = 0;i<lights.size();i++){
        Vec3f light_dir = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        Vec3f shadow_orig = light_dir*N < 0 ? point - N*5e-2 : point + N*5e-2; // checking if the point lies in the shadow of the lights[i]
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        Vec3f tmp_vec;
        int tmp_mat = -1;
        if (scene_intersect(shadow_orig, light_dir, shadow_pt, shadow_N, tmpmaterial,-1,tmp_mat,tmp_vec,-1) && (shadow_pt-shadow_orig).norm() < light_distance){
            if (tmp_mat==4){
                diffuse_light_intensity  += (lights[i].intensity * std::max(0.4f, light_dir*N))*0.4;
                specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
            }
            continue;
        }
        diffuse_light_intensity  += lights[i].intensity * std::max(0.4f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;


    }
    //расчет доп света от огня
    Vec3f fin_fire_color;
    {
        Vec3f fire_color(1,1,1);
        Vec3f fire_center(0,0,-3);
        Vec3f fire_dir = (fire_center-point).normalize();
        Vec3f new_point;
        Vec3f fire_N;
        Material fire_mat;
        int tmp_mat=4;
        scene_intersect(fire_center, fire_dir, new_point, fire_N, fire_mat, 0 ,tmp_mat,fire_color);
        float dist = pow((fire_center-point).norm(),2);
        fire_color =fire_color*((float) 8/dist);
        //Vec3f res = Vec3f(1.1, 1, 1)*1.2;
        //тени
        Vec3f light_dir = (fire_center - point).normalize();
        float light_distance = (-fire_center + point).norm();

        Vec3f shadow_orig = light_dir*N < 0 ? point - N*5e-2 : point + N*5e-2; // checking if the point lies in the shadow of the lights[i]
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        Vec3f tmp_vec;
        tmp_mat = -1;
        if (scene_intersect(shadow_orig, light_dir, shadow_pt, shadow_N, tmpmaterial,-1,tmp_mat,tmp_vec,-1) ){
            if (tmp_mat==4){
                fin_fire_color = Vec3f (material.diffuse_color.x*fire_color.x,material.diffuse_color.y*fire_color.z,material.diffuse_color.x*fire_color.z);
                fin_fire_color = fin_fire_color*material.albedo[0];
            }else{
                fin_fire_color = Vec3f (0,0,0);
            }
        }else{
            fin_fire_color = Vec3f (0,0,0);
        }

    }
    fin_fire_color =fin_fire_color+Vec3f (0.07,0.03,0.01);
    if (cur_mat == -1 && out_m == 4){
        float r_0 = pow(abs((float)(n1-n2)/(n1+n2)),2);
        float cosALPHA = abs(N*dir);
        float r = r_0+(1-r_0)*pow(1-cosALPHA,5);
        return  material.diffuse_color * diffuse_light_intensity * material.albedo[0]+
                Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1]+
                reflect_color*r+
                refract_color*(0.6-r)+
                fin_fire_color;
    }
    if (cur_mat == 4 && out_m == -1){
        float r_0 = pow(abs((float)(n1-n2)/(n1+n2)),2);
        float cosALPHA = abs(N*dir);
        float r = r_0+(1-r_0)*pow(1-cosALPHA,5);
        return  material.diffuse_color * diffuse_light_intensity * material.albedo[0]+
                Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1]+
                reflect_color*r+
                refract_color*(0.6-r)+
                fin_fire_color;
    }
    return  material.diffuse_color * diffuse_light_intensity * material.albedo[0]+
            Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1]+
            reflect_color*material.albedo[2]+
            refract_color*material.albedo[3]
            +fin_fire_color;
}
void render() {
    const int width    = 1024;
    const int height   = 1024;
    const float fov      = M_PI/2.;
    std::vector<Vec3f> framebuffer(width*height);


#pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);//
            Vec3f dir = Vec3f(x, y, -1).normalize(); // вектор запуска луча из камеры
            framebuffer[i+j*width] = ray_marching(Vec3f(0, 0, 0), dir, -1);
        }
        std::cout << j << std::endl;
    }
    std::ofstream ofs; // save the framebuffer to file
    ofs.open("../out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);

        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, c[j])));
        }
    }
    ofs.close();
    std::vector<unsigned char> pixmap(width*height*3);
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        pixmap[i*3] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i].x)));
        pixmap[i*3+1] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i].y)));
        pixmap[i*3+2] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i].z)));

    }
    stbi_write_jpg("../out.jpg", width, height, 3, pixmap.data(), 100);
}
int main() {
    int n = -1;
    unsigned char *pixmap = stbi_load("../map_1.jpg", &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
        return -1;
    }
    envmap = std::vector<Vec3f>(envmap_width*envmap_height);
    for (int j = envmap_height-1; j>=0 ; j--) {
        for (int i = 0; i<envmap_width; i++) {
            envmap[i+j*envmap_width] = Vec3f(pixmap[(i+j*envmap_width)*3+0], pixmap[(i+j*envmap_width)*3+1], pixmap[(i+j*envmap_width)*3+2])*(1/255.);
        }
    }
    stbi_image_free(pixmap);
    render();
    return 0;
}