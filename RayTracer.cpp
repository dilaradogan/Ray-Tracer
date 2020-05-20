#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
using namespace std;

const float EPS = 0.0001;
const float M_EPS = -0.0001;
const int MAX = 4;
const float T_FIRST = 1000.0;
typedef struct {
	float pigNum;
	float surfNum;
	float x;
	float y;
	float z;
	float r;
} Sphere;

typedef struct {
	float r;
	float g;
	float b;
} Pigment;

typedef struct {
	float ka;
	float kd;
	float ks;
	float alpha;
	float kr;
	float kt;
	float n2;
} SurfaceFinish;

typedef struct {
	float x;
	float y;
	float z;
	float r;
	float g;
	float b;
	float ak;
	float bk;
	float ck;
} Light;

typedef struct {
	float x;
	float y;
	float z;
} Point;

typedef struct {
	Point O;
	Point d;
	float t;
} Ray;

typedef struct {
	float i;
	float j;
	float r;
	float g;
	float b;
} Pixel;

int width = 0;
int height = 0;
string output_file;
Point camera;
Point at;
Point up;
int fovy;

int num_lights = 0;
vector<Light> lights;
int numP = 0;
vector<Pigment> pigments;
int numF = 0;
vector<SurfaceFinish> finishes;
int num_sphere = 0;
vector <Sphere> spheres;

Point surface_normal;
int closest_sphere;
ofstream outfile;
bool refraction;

Point cross(Point a, Point b) {
	return Point{ a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x };
}
float clamp(float n, float lower, float upper) {
	return std::max(lower, std::min(n, upper));
}
float dot(Point u, Point v) {
	return u.x * v.x + u.y * v.y + u.z * v.z;
}
float length(Point v) {
	return sqrt(dot(v, v));
}
Point normalize(Point v) {
	return { v.x / length(v), v.y / length(v), v.z / length(v) };
}

Point pointAdd(Point p1, Point p2) {
	return { p1.x + p2.x, p1.y + p2.y, p1.z + p2.z };
}
Point pointMinus(Point p1, Point p2) {
	return { p1.x - p2.x, p1.y - p2.y, p1.z - p2.z };
}
Point pointMultiply(Point p1, Point p2) {
	return { p1.x * p2.x, p1.y * p2.y, p1.z * p2.z };
}
Point scalerProduct(Point p, float f) {
	return { p.x * f, p.y * f, p.z * f };
}
Ray compute(int i, int j) {
	float aspect_ratio = (width * 1.0) / height;
	Point co = camera;
	Point cz = scalerProduct(normalize(pointMinus(at, camera)), -1.0);
	Point cx = normalize(cross(up, cz));
	Point cy = cross(cz, cx);
	float h = 2 * tan(fovy / 2.0);//2 * tan(fovy / 2);
	float w = h * aspect_ratio;
	Point p;
	p.x = ((w * ((2 * j) - width)) / (2.7 * width));
	p.y = ((h * ((-2 * i) + height)) / (2.7 * height));
	p.z = -1.0;
	Point X = scalerProduct(cx, p.x);
	Point Y = scalerProduct(cy, p.y);
	Point Z = scalerProduct(cz, p.z);
	p = pointAdd(pointAdd(co, X), pointAdd(Y, Z));
	Ray ray;
	ray.O = co;
	ray.d = normalize(pointMinus(p, co));
	ray.t = T_FIRST;
	/*cout << "compute ray R(i,j): " << "\n\tr.O: " << ray.O.x << ", " << ray.O.y << ", " << ray.O.z << "\n";
	cout << "\tr.d: " << ray.d.x << ", " << ray.d.y << ", " << ray.d.z << "\n";
	cout << "\tr.t: " << ray.t<< "\n";*/
	return ray;
}

Ray lightObjectIntersection(Sphere sphere, Ray ray) {
	float r = sphere.r;
	Point c = { sphere.x,sphere.y,sphere.z };
	Point o = ray.O;
	Point u = { c.x - o.x,c.y - o.y,c.z - o.z };
	Point d = ray.d;
	float d_dot_u = dot(d, u);
	Point x = pointMinus(u, scalerProduct(d, d_dot_u));
	float delta = 4 * ((r * r) - (length(x) * length(x)));
	float t1, t2;
	if (delta < EPS) {
		return ray;
	}
	else {
		t1 = d_dot_u - sqrt((r * r) - (length(x) * length(x)));
		t2 = d_dot_u + sqrt((r * r) - (length(x) * length(x)));
	}
	if (t1 > EPS&& t2 > EPS) {
		if (t1 > t2) {
			if (ray.t > t2) {
				ray.t = t2;
			}
		}
		else {
			if (ray.t > t1) {
				ray.t = t1;
			}
		}
	}
	else if (t1 > EPS&& t2 < EPS) {
		if (ray.t > t1) {
			ray.t = t1;
		}
	}
	else if (t2 > EPS&& t1 < EPS) {
		if (ray.t > t2) {
			ray.t = t2;
		}
	}
	return ray;
}
Point pointDiv(Point p1, Point p2) {
	return { p1.x / p2.x, p1.y / p2.y, p1.z / p2.z };
}
int visible(Point p, int light_source) {
	Ray r;
	r.O = { lights[light_source].x, lights[light_source].y, lights[light_source].z };
	r.d = normalize(pointMinus(p, r.O));
	r.t = T_FIRST;
	Ray ray = r;
	for (int i = 0; i < num_sphere; i++) {
		ray = lightObjectIntersection(spheres[i], ray);
	}
	Point pNew = pointAdd(ray.O, scalerProduct(ray.d, ray.t));
	Point pa = pointMinus(p, pNew);
	Point pb = pointMinus(pNew, p);
	if ((pa.x<EPS && pa.x>M_EPS&&pa.y<EPS &&pa.y>M_EPS && pa.z<EPS && pa.z>M_EPS)|| (pb.x<EPS && pb.x>M_EPS && pb.y<EPS && pb.y>M_EPS && pb.z<EPS && pb.z>M_EPS)){
		return 1;
	}
	else
		return 0;
}

Ray raySphereIntersection(Sphere sphere, Ray& ray, int sphere_num) {
	float r = sphere.r;
	Point c = { sphere.x,sphere.y,sphere.z };
	Point o = ray.O;
	Point u = { c.x - o.x,c.y - o.y,c.z - o.z };
	Point d = ray.d;
	float d_dot_u = dot(d, u);
	Point x = pointMinus(u, scalerProduct(d, d_dot_u));
	float delta = 4 * ((r * r) - (length(x) * length(x)));
	float t1, t2;
	if (delta < EPS) {
		return ray;
	}
	else {
		t1 = d_dot_u - sqrt((r * r) - (length(x) * length(x)));
		t2 = d_dot_u + sqrt((r * r) - (length(x) * length(x)));
	}
	if (t1 > EPS&& t2 > EPS) {
		if (t1 > t2) {
			if (ray.t > t2) {
				ray.t = t2;
				surface_normal = normalize(pointMinus(pointAdd(ray.O, scalerProduct(ray.d, ray.t)), c));
				closest_sphere = sphere_num;
			}
		}
		else {
			if (ray.t > t1) {
				ray.t = t1;
				surface_normal = normalize(pointMinus(pointAdd(ray.O, scalerProduct(ray.d, ray.t)), c));
				closest_sphere = sphere_num;
			}
		}
	}
	else if (t1 > EPS&& t2 < EPS) {
		if (ray.t > t1) {
			ray.t = t1;
			surface_normal = normalize(scalerProduct(pointMinus(pointAdd(ray.O, scalerProduct(ray.d, ray.t)), c), -1));
			closest_sphere = sphere_num;
		}
	}
	else if (t2 > EPS&& t1 < EPS) {
		if (ray.t > t2) {
			ray.t = t2;
			surface_normal = normalize(scalerProduct(pointMinus(pointAdd(ray.O, scalerProduct(ray.d, ray.t)), c), -1));
			closest_sphere = sphere_num;
		}
	}
	return ray;
}
Ray intersect(Ray ray) {
	for (int i = 0; i < num_sphere; i++) {
		ray = raySphereIntersection(spheres[i], ray, i);
	}
	return ray;
}
Point phong(int light_num, Point p, Point normal) {
	float di = length(pointMinus(p, { lights[light_num].x, lights[light_num].y, lights[light_num].z }));
	float c_factor = di * di * lights[light_num].ck;
	float b_factor = di * lights[light_num].bk;
	float div_factor = lights[light_num].ak + b_factor + c_factor;

	int pig_num = spheres[closest_sphere].pigNum;
	Point color = { pigments[pig_num].r, pigments[pig_num].g,  pigments[pig_num].b };
	int surf_num = spheres[closest_sphere].surfNum;
	Point l = normalize(pointMinus({ lights[light_num].x,  lights[light_num].y, lights[light_num].z }, p));
	Point v = normalize(pointMinus(camera, p));
	Point h = normalize(pointAdd(l, v));

	float Kd = finishes[surf_num].kd;
	float Kddot = ((dot(normal, l) > EPS) ? dot(normal, l) : EPS);
	float Ks = finishes[surf_num].ks;
	double Ksdot = pow((dot(normal, h) > 0.0) ? dot(normal, h) : 0.0, finishes[surf_num].alpha);
	Point Li = { lights[light_num].r, lights[light_num].g, lights[light_num].b };
	Point Li2 = Li;
	Li = pointMultiply(Li, color);
	Point pkd = scalerProduct(scalerProduct(Li, Kd), Kddot);
	Point pks = scalerProduct(Li2, Ks * Ksdot);
	//cout << "pks: " << pks.x << ", " << pks.y << ", " << pks.z << "\n";
	Point Lii = pointAdd(pkd, pks);
	Lii = scalerProduct(Lii, 1.0 / div_factor);
	Lii.x = clamp(Lii.x, 0.0, 1.0);
	Lii.y = clamp(Lii.y, 0.0, 1.0);
	Lii.z = clamp(Lii.z, 0.0, 1.0);
	//cout << "Lii: " << Lii.x << ", " << Lii.y << ", " << Lii.z << "\n";
	return Lii;
}
Ray transmit(Ray R, Point p, Point normal) {
	Point v = scalerProduct(normalize(R.d), -1.0);
	Ray Rt;
	Rt.O = p;
	float vn = dot(normal, v);
	int surf_num = spheres[closest_sphere].surfNum;
	float n2 = finishes[surf_num].n2;
	float n = 0.29 / n2;
	float x = sqrt(1.0-(pow(n,2)*(1-(pow(vn, 2)))));
	float a = (n * vn) - x;
	Point d = scalerProduct(normal, a);
	d = pointMinus(d, scalerProduct(v, n));
	Rt.d = d;
	Rt.t = T_FIRST;
	return Rt;
}
Ray reflect(Ray R, Point p, Point normal) {
	Point v = scalerProduct(normalize(R.d), -1.0);
	Ray Rr;
	Rr.O = p;
	float nv = 2 * dot(normal, v);
	Point partial = scalerProduct(normal, nv);
	partial = pointMinus(partial, v);
	Rr.d = partial;
	Rr.t = T_FIRST;
	return Rr; 
}
Point trace(Ray R, int depth) {
	Point localC, reflectedC, transmittedC;
	Point p;
	Point normal;
	if (depth > MAX) {
		return { 0.5, 0.5, 0.5 };
	}
	R = intersect(R);
	p = pointAdd(R.O, scalerProduct(R.d, R.t));
	if (R.t == T_FIRST) {
		return { 0.5, 0.5, 0.5 };
	}
	normal = surface_normal;
	localC = { 0,0,0 };
	int pig_num = spheres[closest_sphere].pigNum;
	int surf_num = spheres[closest_sphere].surfNum;
	Point colorC = { pigments[pig_num].r, pigments[pig_num].g,  pigments[pig_num].b };
	//cout << colorC.x<< ", "<<colorC.y <<", "<< colorC.z << "\n";
	Point La = { lights[0].r, lights[0].g, lights[0].b };
	localC = pointMultiply(scalerProduct(La, finishes[surf_num].ka),colorC);
	for (int i = 1; i < num_lights; i++) {
		if (visible(p, i) == 1) {
			localC = pointAdd(localC, phong(i, p, normal));
		}
	}

	reflectedC = { 0,0,0 };
	float kr = finishes[surf_num].kr;
	if (kr > EPS) {//Reflective surface
		Ray Rr = reflect(R, p, normal);
		reflectedC = trace(Rr, depth + 1);
	}
	transmittedC = { 0,0,0 };
	float kt = finishes[surf_num].kt;
	if (kt > EPS) {//Transparent surface
		Ray Rt = transmit(R, p, normal);
		transmittedC = trace(Rt, depth + 1);
	}

	localC = pointAdd(localC, scalerProduct(reflectedC, kr));
	localC = pointAdd(localC, scalerProduct(transmittedC, kt));
	localC.x = clamp(localC.x, 0.0, 1.0);
	localC.y = clamp(localC.y, 0.0, 1.0);
	localC.z = clamp(localC.z, 0.0, 1.0);
	
	return localC;
}

void write_pixel(int i, int j, Point color) {
	float r = color.x;
	float g = color.y;
	float b = color.z;
	int ir = int(255.99 * r);
	int ig = int(255.99 * g);
	int ib = int(255.99 * b);
	outfile << ir << " " << ig << " " << ib << "\n";
}
int
main(int argc, char** argv)
{
	int refl = 1;
	if (argc > 2) {
		refl = strcmp(argv[2], "y");
	}
	if (refl== 0) {
		refraction = 1;
	}
	else {
		refraction = 0;
	}
	width = 0;
	height = 0;
	char s[20];
	FILE* fp;
	vector <string> tokens;
	fp = fopen(argv[1], "r");
	char out_file[100] = "";
	if (fp == NULL) {
		printf("No file found.");
		return 0;
	}
	int i = 0;
	fscanf(fp, "%s", &out_file);
	while (fscanf(fp, "%s", s) == 1) {
		tokens.push_back(s);
	}
	int index = 0;
	width = stoi(tokens[index++]);
	height = stoi(tokens[index++]);
	camera = { stof(tokens[index++]),stof(tokens[index++]),stof(tokens[index++]) };
	at = { stof(tokens[index++]),stof(tokens[index++]),stof(tokens[index++]) };
	up = { stof(tokens[index++]),stof(tokens[index++]),stof(tokens[index++]) };
	fovy = stoi(tokens[index++]);
	num_lights = stoi(tokens[index++]);
	for (int i = 0; i < num_lights; i++) {
		lights.push_back({ stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]) });
	}
	numP = stoi(tokens[index++]);
	for (int i = 0; i < numP; i++) {
		index++;
		pigments.push_back({ stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]) });
	}
	numF = stoi(tokens[index++]);
	for (int i = 0; i < numF; i++) {
		if (refraction == true) {
			finishes.push_back({ stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]) });
		}
		else {
			finishes.push_back({ stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), stof(tokens[index++]), 0.0 , 0.0 });
		}
	}
	num_sphere = stoi(tokens[index++]);
	for (int i = 0; i < num_sphere; i++) {
		Sphere sp;
		sp.pigNum = stof(tokens[index++]);
		sp.surfNum = stof(tokens[index++]);
		index++;
		sp.x = stof(tokens[index++]);
		sp.y = stof(tokens[index++]);
		sp.z = stof(tokens[index++]);
		sp.r = stof(tokens[index++]);
		spheres.push_back(sp);
	}
	cout << out_file << "\n" << width << "	" << height << "\n";
	cout << camera.x << "	" << camera.y << "	" << camera.z << "\n";
	cout << at.x << "	" << at.y << "	" << at.z << "\n";
	cout << up.x << "	" << up.y << "	" << up.z << "\n";
	cout << fovy << "\n" << num_lights << "\n";
	for (int i = 0; i < num_lights; i++) {
		cout << lights[i].x << " " << lights[i].y << " " << lights[i].z << " " << lights[i].r << " " << lights[i].g << " " << lights[i].b << " " << lights[i].ak << " " << lights[i].bk << " " << lights[i].ck << "\n";
	}
	cout << numP << "\n";
	for (int i = 0; i < numP; i++) {
		cout << pigments[i].r << " " << pigments[i].g << " " << pigments[i].b << "\n";
	}
	cout << numF << "\n";
	for (int i = 0; i < numF; i++) {
		cout << finishes[i].ka << " " << finishes[i].kd << " " << finishes[i].ks << " " << finishes[i].alpha << " " << finishes[i].kr << " " << finishes[i].kt << " " << finishes[i].n2 <<"\n";
	}
	cout << num_sphere << "\n";
	for (int i = 0; i < num_sphere; i++) {
		cout << spheres[i].pigNum << " " << spheres[i].surfNum << " " << spheres[i].x << " " << spheres[i].y << " " << spheres[i].z << " " << spheres[i].r << "\n";
	}
	fclose(fp);

	outfile.open(out_file, ios_base::out);
	outfile << "P3\n" << width << " " << height << "\n255\n";
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			Ray R = compute(i, j);
			Point color = trace(R, 0);
			write_pixel(i, j, color);
		}
	}
	outfile.close();

	return 0;
}