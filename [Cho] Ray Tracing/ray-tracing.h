#include <list>
using namespace std;

class Color;
class Light;
class Figure;
class Plane;
class Sphere;
class Ray;

class Vec
{
  protected:
    double x;
    double y;
    double z;

  public:
    Vec();
    Vec(ifstream& ifs);
	Vec(double x, double y, double z);
	const double dotProduct(Vec& thing) const;
	const Vec subtract(Vec otherVec) const;
	const Vec add(Vec otherVec) const;
	const Vec scale(double t) const;
	const Vec divide(double d) const;
	const double getLength() const;
	const Vec makeNormal() const;
	
};

class Ray
{
protected:
	Vec start;
	Vec end;

public:
	Ray();
	Ray(Vec& p0, Vec& p1);
	const Vec getStart() const;
	const Vec getEnd() const;
};


class Color
{
  friend Color operator*(double num, const Color& c);

  protected:
   double red, green, blue;

  public:
	Color();
    Color(ifstream& ifs);
    Color(double r, double g, double b);
	double getRed();
	double getGreen();
	double getBlue();
	Color add(Color otherColor);
	Color multiply(Color otherColor);
	Color scale(double t);
	Color mini();
	bool isEqualTo(Color otherColor);
};

class Light
{
  public:
    Light(ifstream& ifs);
	const Vec getPosition() const;
	Color getShading();
  private:
    Vec position;
    Color shading;
    double c0, c1, c2;
};


class Figure
{
 protected:
   Color ambient;
   Color diffuse;
   Color specular;
   Color reflectivity;
   Color transmissivity;
   double shininess;
   double indexOfRefraction;
   int rFlag, tFlag;

 public:
   Figure();
   void initFigure(ifstream& ifs);
   virtual double intersection(const Ray& r, double minT, double maxT) const = 0;
   virtual Vec* getNormal(Vec* p) = 0;
   Color getAmbient();
   Color getDiffuse();
   Color getSpecular();
   Color getReflectivity();
   Color getTransmissivity();
   double getShininess();
   double getRefraction();
};



class Plane : public Figure
{
  private:
    Vec abcVector;
    double dScalar;
    Vec direction1;
    Vec direction2;
  public:
    Plane(ifstream& ifs);
	double intersection(const Ray& r, double minT, double maxT) const;
	Vec* getNormal(Vec* p);
};

class Sphere : public Figure
{
  private:
    Vec center;
    double radius;
  public:
    Sphere(ifstream& ifs);
	double intersection(const Ray& r, double minT, double maxT) const;
	Vec* getNormal(Vec* p);
};

