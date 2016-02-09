#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <cassert>
#include <list>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

#include "ray-tracing.h"

//vector<vector<Color> > image;

Color RT_trace(const Ray& r, double depth);
Color RT_shade(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal, bool entering, double depth);
Color RT_transmit(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal, bool entering, double depth);
Color RT_reflect(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal, double depth);

list<Figure*> shapeList;
list<Light*> lightList;

Color image[1000][1000];

double cameraX, cameraY, cameraZ;
int horizontalResolution, verticalResolution;
double zCoor;
double minX, maxX;
double minY, maxY;

int maxDepth;
Color backgroundColor;
Color ambient;

Color pixColor;
const double MinT = .01;
const double MaxT = 1000000;

Vec::Vec() : x(0.0),y(0.0),z(0.0) {}

Vec::Vec(ifstream& ifs)
{
  ifs >> x >> y >> z;
}

Vec::Vec(double x, double y, double z) : x(x), y(y), z(z)
{}

const double Vec::dotProduct(Vec& thing) const
{
	double dpx, dpy, dpz, total;
	dpx = x * thing.x;
	dpy = y * thing.y;
	dpz = z * thing.z;
	total = dpx + dpy + dpz;
	return total;
}

const Vec Vec::subtract(Vec otherVec) const
{
	Vec result = Vec();
	result.x = x - otherVec.x;
	result.y = y - otherVec.y;
	result.z = z - otherVec.z;
	return result;
}

const Vec Vec::add(Vec otherVec) const
{
	Vec result = Vec();
	result.x = x + otherVec.x;
	result.y = y + otherVec.y;
	result.z = z + otherVec.z;
	return result;
}

const Vec Vec::scale(double t) const
{
	Vec result = Vec();
	result.x = x * t;
	result.y = y * t;
	result.z = z * t;
	return result;
}

const Vec Vec::divide(double d) const
{
	Vec result = Vec();
	result.x = x / d;
	result.y = y / d;
	result.z = z / d;
	return result;
}

const double Vec::getLength() const
{
	return sqrt((x * x) + (y * y) + (z * z));
}

const Vec Vec::makeNormal() const
{
	return this->divide(this->getLength());
}

Ray::Ray() : start(Vec()), end(Vec())
{
}

Ray::Ray(Vec& p0, Vec& p1) : start(p0), end(p1)
{
}

const Vec Ray::getStart() const
{
	return start;
}

const Vec Ray::getEnd() const
{
	return end;
}

Color::Color(): red(0.0), green(0.0), blue(0.0) {};

Color::Color(ifstream& ifs)
{
  ifs >> red >> green >> blue;
}

Color::Color(double r, double g, double b) : red(r), green(g), blue(b)
{}

double Color::getRed()
{
	return red;
}

double Color::getGreen()
{
	return green;
}

double Color::getBlue()
{
	return blue;
}

Color Color::add(Color otherColor)
{
	Color result = Color();
	result.red = red + otherColor.red;
	result.green = green + otherColor.green;
	result.blue = blue + otherColor.blue;
	return result;
}

Color Color::multiply(Color otherColor)
{
	Color result = Color();
	result.red = red * otherColor.red;
	result.green = green * otherColor.green;
	result.blue = blue * otherColor.blue;
	return result;
}

Color Color::scale(double t)
{
	Color result = Color();
	result.red = red * t;
	result.green = green * t;
	result.blue = blue * t;
	return result;
}

Color Color::mini()
{
	Color result = Color();
	result.red = min(red, 1.0);
	result.green = min(green, 1.0);
	result.blue = min(blue, 1.0);
	return result;
}

bool Color::isEqualTo(Color otherColor)
{
	if (red == otherColor.red && green == otherColor.green && blue == otherColor.blue)
	{
		return true;
	}
	else return false;
}

Figure::Figure(){}

void Figure::initFigure(ifstream& ifs)
{
 ambient = Color(ifs);
 diffuse = Color(ifs);
 specular = Color(ifs);
 reflectivity = Color(ifs);
 transmissivity = Color(ifs);
 ifs >> shininess >> indexOfRefraction >> rFlag >> tFlag;
}

Color Figure::getAmbient()
{
	return ambient;
}

Color Figure::getDiffuse()
{
	return diffuse;
}

Color Figure::getSpecular()
{
	return specular;
}

Color Figure::getReflectivity()
{
	return reflectivity;
}

Color Figure::getTransmissivity()
{
	return transmissivity;
}

double Figure::getShininess()
{
	return shininess;
}

double Figure::getRefraction()
{
	return indexOfRefraction;
}

Light::Light(ifstream& ifs) : position(ifs), shading(ifs)
{
  ifs >> c0 >> c1 >> c2;
}

const Vec Light::getPosition() const
{
	return position;
}

Color Light::getShading()
{
	return shading;
}

Sphere::Sphere(ifstream& ifs) : center(ifs)
{
  ifs >> radius;
  initFigure(ifs);
}

Plane::Plane(ifstream& ifs) : abcVector(ifs)
{
  ifs >> dScalar;
  initFigure(ifs);
  direction1 = Vec(ifs);
  direction2 = Vec(ifs);
}

void parseSceneFile(char* sceneName)
{
  double bgr, bgg, bgb;
  double ar, ag, ab;
  ifstream ifs;
  assert (sceneName != 0);
  ifs.open(sceneName);
  ifs >> cameraX;
  ifs >> cameraY;
  ifs >> cameraZ;
  ifs >> zCoor;
  ifs >> minX >> maxX;
  ifs >> minY >> maxY;
  ifs >> bgr >> bgg >> bgb;
  backgroundColor = Color(bgr,bgg,bgb);
  ifs >> ar >> ag >> ab;
  ambient = Color(ar,ag,ab);
  ifs >> maxDepth;
  ifs >> horizontalResolution >> verticalResolution;
  int numLights, numSpheres, numPlanes;
  ifs >> numLights;
  for (int i=0; i<numLights; ++i) lightList.push_front(new Light(ifs));
  ifs >> numSpheres;
  for (int i=0; i<numSpheres; ++i) shapeList.push_front(new Sphere(ifs));
  ifs >> numPlanes;
  for (int i=0; i<numPlanes; ++i) shapeList.push_front(new Plane(ifs));
  ifs.close();
}

Vec* Sphere::getNormal(Vec* p)
{
	return new Vec(p->subtract(center));
}

Vec* Plane::getNormal(Vec* p)
{
	return new Vec(abcVector.makeNormal());
}

double Sphere::intersection(const Ray& r, double minT, double maxT) const
{
	Vec h = r.getEnd().subtract(r.getStart());
	Vec k = r.getStart().subtract(center);

	double a = h.dotProduct(h);
	double b = 2 * h.dotProduct(k);
	double c = k.dotProduct(k) - (radius * radius);

	double discri = (b*b) - (4 * a * c);
	if (discri > 0) //two roots
	{
		double t1 = (-b + sqrt(discri)) / (2 * a);
		double t2 = (-b - sqrt(discri)) / (2 * a);
		if (t1 > t2 && t2 > 0)
		{
			return t2;
		}
		else if (t1 < t2 && t1 > 0) return t1;
	}
	else if (discri == 0) //one root
	{
		double t = -b / (2 * a);
		return t;
	}
	else //if (discri < 0) //no intersection
	{
		return maxT + 2;
	}
}

double Plane::intersection(const Ray& r, double minT, double maxT) const
{
	//D = p dotproduct n
	//t = (D - r.getStart() dotproduct n)/((r.getEnd() - r.getStart()) dotproduct n) //:(

	//Ray rr = r;
	Vec p10 = (r.getEnd()).subtract(r.getStart());
	Vec n = abcVector;
	//Vec p = r.getStart(); //?
	double D = n.dotProduct(p10);
	double deno = p10.dotProduct(n);
	double t = (D - ((r.getStart()).dotProduct(n))) / deno;
	if (deno == 0) //parallel to plane
	{
		return maxT + 2;
	}
	else //if (deno != 0) //hits plane
	{
		return t;
	}
}

pair<double, Figure*> nearestIntersection(const Ray& r,
	double minT, double maxT,
	bool mayBeTransparent = true)
{
	pair<double, Figure*> blah;
	double thing = 0.0;
	Figure* figgy = NULL;
	for (list<Figure*>::const_iterator iter = shapeList.cbegin();
	iter != shapeList.cend();
		++iter)
	{
		double smallest = maxT + 1;
		double inter = (*iter)->intersection(r, minT, maxT);
		if ((inter < thing || thing <= 0.0) && inter > minT && inter < smallest)
		{
			thing = inter;
			figgy = *iter;
			blah = make_pair(thing, figgy);
		}
	}
	return blah;
}

Color specularShade(Figure* obj, const Vec& normal,
	Light* light, const Vec& lightDirection, double dotProduct,
	const Ray& ray)
{
	//math?
	Color i = light->getShading();
	Color k = obj->getSpecular();
	double sr = i.getRed() * k.getRed();
	double sg = i.getGreen() * k.getGreen();
	double sb = i.getBlue() * k.getBlue();
	double p = max(dotProduct, 0.0);
	double po = pow(p, obj->getShininess());
	Color newColor = Color(sr, sg, sb);
	Color newNewColor = newColor.scale(po);
	return newNewColor;
}

Color diffuseShade(Figure* obj, Light* light, double dotProduct)
{
	//math?
	Color i = light->getShading();
	Color k = obj->getDiffuse();
	double dr = i.getRed() * k.getRed();
	double dg = i.getGreen() * k.getGreen();
	double db = i.getBlue() * k.getBlue();
	double p = max(dotProduct, 0.0);
	Color newColor = Color(dr, dg, db);
	Color newNewColor = newColor.scale(p);
	return newNewColor;
}

Color RT_lights(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal)
{
	Color newColor = Color();
	for (list<Light*>::const_iterator iter = lightList.cbegin();
	iter != lightList.cend();
		++iter)
	{
		Light* s = *iter;
		Vec ii = i;
		Vec ss = s->getPosition();
		Vec nn = normal;
		Ray L_ray = Ray(ii, ss);
		Vec L = s->getPosition().subtract(i).makeNormal();
		pair<double, Figure*> inter = nearestIntersection(L_ray, MinT, 1.0, false);
		if (inter.first == 0)
		{
			double dotPro = max(L.dotProduct(nn), 0.0);
			newColor = newColor.add(diffuseShade(obj, s, dotPro));
			newColor = newColor.add(specularShade(obj, nn, s, L, dotPro, ray));
		}
	}
	return newColor;
}

Color RT_shade(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal, bool entering, double depth)
{
	Color newColor = ambient.multiply(obj->getAmbient());
	newColor = newColor.add(RT_lights(obj, ray, i, normal));
	if (depth < maxDepth)
	{
		newColor = newColor.add(RT_reflect(obj, ray, i, normal, depth));
		newColor = newColor.add(RT_transmit(obj, ray, i, normal, entering, depth));
	}
	return newColor;

}

Color RT_trace(const Ray& r, double depth)
{
	pair<double, Figure*> i = nearestIntersection(r, MinT, MaxT);
	double t = i.first;
	Figure* figure = i.second;

	//for testing first two images
	//if (figure == NULL) {
	//	return backgroundColor;
	//}
	//else
	//{
	//	return figure->getAmbient();
	//}

	if (figure == NULL)
	{
		return backgroundColor;
	}
	else
	{
		Vec p = (r.getStart()).add((r.getEnd().subtract(r.getStart())).scale(t));
		Vec* n = figure->getNormal(&p);
		Vec norm = n->makeNormal();
		return RT_shade(figure, r, p, norm, false, depth).mini();
	}
}

Color RT_transmit(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal,
	bool entering, double depth)
{
	Color nothing = Color();
	if (!(obj->getTransmissivity().isEqualTo(nothing))) //does have transmissive qualities
	{
		Ray rr = ray;
		Vec nn = normal;
		Vec inter = i;

		Vec q = nn.scale(cos(obj->getRefraction()));
		Vec s = q.subtract(i);
		Vec m = s.divide(sin(obj->getRefraction()));
		Vec p = (nn.scale(-1)).scale(cos(obj->getRefraction()));
		Vec w = m.scale(obj->getRefraction());
		Vec refracting = p.add(w);
		Ray T_ray = Ray(inter, refracting);
		Vec T = refracting;
		if (T.dotProduct(nn) < 0)
		{
			Color newColor = RT_trace(T_ray, depth + 1);
			newColor = newColor.multiply(obj->getTransmissivity());
			return newColor;
		}
		else return Color();
	}
	else return Color();
}

Color RT_reflect(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal, double depth)
{
	Color nothing = Color();
	if (!(obj->getReflectivity().isEqualTo(nothing))) //does have reflective qualities
	{
		Vec nn = normal;
		Vec inter = i;
		Vec v = i.scale(-1);
		Vec q = nn.scale(nn.dotProduct(v));
		Vec s = q.subtract(v);
		Vec r = v.add(s.scale(2));
		Ray R_ray = Ray(inter, r);
		Color newColor = RT_trace(R_ray, depth + 1);
		newColor = newColor.multiply(obj->getReflectivity());
		return newColor;
	}
	else return Color();
}

void RT_algorithm()
{
	for (int j = 0; j < verticalResolution; ++j)
	{
		for (int i = 0; i < horizontalResolution; ++i)
		{
			double w = (maxX - minX) / horizontalResolution;
			double h = (maxY - minY) / verticalResolution;
			double x = minX + i * w + w / 2;
			double y = minY + j * h + h / 2;
			Vec p11 = Vec(x, y, zCoor); //zCoor will usually be 0 for this assignment
			Vec p00 = Vec(cameraX, cameraY, cameraZ);
			Ray R = Ray(p00, p11);
			pixColor = RT_trace(R, 1);
			image[i][verticalResolution - j] = pixColor;
		}
	}
}

void writeImageFile()
{
	ofstream ofs("picture.ppm");
	ofs << "P3" << endl;
	ofs << horizontalResolution << " " <<  verticalResolution << endl;
	ofs << "255" << endl;
	for (int j = 0; j < verticalResolution; ++j)
	{
		for (int i = 0; i < horizontalResolution; ++i)
		{
			double r = image[i][j].getRed();
			double g = image[i][j].getGreen();
			double b = image[i][j].getBlue();
			r = (int)(r * 255);
			g = (int)(g * 255);
			b = (int)(b * 255);
			ofs << r << " " << g << " " << b << endl;
		}
	}
	ofs.close();
}

int main(int, char *argv[])
{
    parseSceneFile(argv[1]);
	//initializeImage();
	RT_algorithm();
	writeImageFile();
    return 0;
}

