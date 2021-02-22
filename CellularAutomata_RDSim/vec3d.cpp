/***************************************************************************
                          vec3d.cpp  -  description
                             -------------------
    begin                : Wed May 22 2002
    copyright            : (C) 2002 by Malte Schmick
    email                : schmick@mpi-dortmund.mpg.de
 ***************************************************************************/

#include "vec3d.h"
using namespace std;

ostream& operator<<(ostream& strm, const vec3d& V)
{
	V.printStr(strm);
	return strm;
}

istream& operator>>(istream& strm, vec3d& V)
{
	V.InStr(strm);
	return strm;
}

void vec3d::printStr(ostream& strm) const
{
	strm << x << " " << y << " " << z;
}

void vec3d::InStr(istream& strm)
{
	strm >> x >> y >> z;
	if (! strm) return;
}

double vec3d::X() const
{
	return x;
}

double vec3d::Y() const
{
	return y;
}

double vec3d::Z() const
{
	return z;
}

vec3d::vec3d( double X, double Y, double Z)
{
	x=X;
	y=Y;
	z=Z;
}

vec3d::vec3d( vec3d *V)
{
	x=V->x;
	y=V->y;
	z=V->z;
}

vec3d::vec3d( const vec3d &V)
{
	x=V.x;
	y=V.y;
	z=V.z;
}

vec3d::~vec3d( void)
{
}


vec3d operator*(const double& r,  const vec3d& V)
{
	return vec3d (r*V.x, r*V.y, r*V.z);
}

vec3d vec3d::operator*(const double& r) const
{
	return vec3d (r*x,r*y,r*z);
}

vec3d vec3d::operator/(const vec3d& V) const
{
	return vec3d(y*V.z-z*V.y,
		        z*V.x-x*V.z,
		        x*V.y-y*V.x);
}

vec3d vec3d::operator/(const double& r) const
{
	return vec3d (x/r,y/r,z/r);
}

vec3d vec3d::operator+(const vec3d& V) const
{
	return vec3d (x+V.x,y+V.y,z+V.z);
}


vec3d vec3d::operator-(const vec3d& V) const
{
	return vec3d (x-V.x,y-V.y,z-V.z);
}

double vec3d::operator*(const vec3d& V) const
{
	return x*V.x+y*V.y+z*V.z;
}

const double vec3d::SQR(void) const
{
	return (x*x+y*y+z*z);
}

double vabs(const vec3d& V)
{
	return sqrt(V.SQR());
}

double sqr(const vec3d& V)
{
	return V.SQR();
}

const int vec3d::operator==(const vec3d& V) const
{
	if (V.x==0.0 && V.y==0.0 && V.z==0.0)  return (sqr(*this)<=grenze);
	else if (x==0.0 && y==0.0 && z==0.0) return (sqr(V)<=grenze);
	return (sqr(*this-V)/sqr(*this)<=grenze);
}

const int vec3d::operator!=(const vec3d& V) const
{
	return (!(*this==V));
}

vec3d vec3d::operator-()
{
	return vec3d (-x,-y,-z);
}

const vec3d vec3d::unit_vector () const
{
	if (sqr(*this)<=grenze) return vec3d(0.0,0.0,0.0);

	else return (*this/vabs(*this));
}

void vec3d::normalize()
{
	double f=sqr(*this);
	if (f<=grenze)
	{
	  x=0.;
	  y=0.;
	  z=0.;
	}
	else *this/=sqrt(f);
}

const vec3d& vec3d::operator+=(const vec3d& V)
{
	x+= V.x;
	y+= V.y;
	z+= V.z;
	return *this;
}

const vec3d& vec3d::operator-=(const vec3d& V)
{
	x-= V.x;
	y-= V.y;
	z-= V.z;
	return *this;
}

const vec3d& vec3d::operator*=(const double& a)
{
	x*= a;
	y*= a;
	z*= a;
	return *this;
}

const vec3d& vec3d::operator/=(const double& a)
{
	x/= a;
	y/= a;
	z/= a;
	return *this;
}
