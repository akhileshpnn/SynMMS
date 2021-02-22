#include "vec2d.h"

ostream& operator<<(ostream& strm, const vec2d& V)
{
	V.printStr(strm);
	return strm;
}

istream& operator>>(istream& strm, vec2d& V)
{
	V.InStr(strm);
	return strm;
}

void vec2d::printStr(ostream& strm) const
{
	strm << x << " " << y << " ";
}

void vec2d::InStr(istream& strm)
{
	strm >> x >> y;
	if (! strm) return;
}

double vec2d::X() const
{
	return x;
}

double vec2d::Y() const
{
	return y;
}

vec2d::vec2d( double X, double Y)
{
	x=X;
	y=Y;
}

vec2d::vec2d( vec2d *V)
{
	x=V->x;
	y=V->y;
}

vec2d::vec2d( const vec2d &V)
{
	x=V.x;
	y=V.y;
}

vec2d::~vec2d( void)
{
}


vec2d operator*(const double& r,  const vec2d& V)
{
	return vec2d (r*V.x, r*V.y);
}

vec2d vec2d::operator*(const double& r) const
{
	return vec2d (r*x,r*y);
}

vec2d vec2d::operator/(const double& r) const
{
	return vec2d (x/r,y/r);
}

vec2d vec2d::operator+(const vec2d& V) const
{
	return vec2d (x+V.x,y+V.y);
}


vec2d vec2d::operator-(const vec2d& V) const
{
	return vec2d (x-V.x,y-V.y);
}

double vec2d::operator*(const vec2d& V) const
{
	return x*V.x+y*V.y;
}

const double vec2d::SQR(void) const
{
	return (x*x+y*y);
}

double vabs(const vec2d& V)
{
	return sqrt(V.SQR());
}

double sqr(const vec2d& V)
{
	return V.SQR();
}

const bool vec2d::operator==(const vec2d& V) const
{
	if (V.x==0.0 && V.y==0.0)  return (sqr(*this)<=grenze);
	else if (x==0.0 && y==0.0) return (sqr(V)<=grenze);
	return (sqr(*this-V)/sqr(*this)<=grenze);
}

const bool vec2d::operator!=(const vec2d& V) const
{
	return (!(*this==V));
}

const bool vec2d::operator  >(const vec2d& V) const
{
	return (sqr(*this)>V.SQR());
}

const bool vec2d::operator  <(const vec2d& V) const
{
	return (sqr(*this)<V.SQR());
}

const bool vec2d::operator >=(const vec2d& V) const
{
	return (sqr(*this)>=V.SQR());
}

const bool vec2d::operator <=(const vec2d& V) const
{
	return (sqr(*this)<=V.SQR());
}

vec2d vec2d::operator-()
{
	return vec2d (-x,-y);
}

vec2d vec2d::operator!()
{
	return vec2d (-y,x);
}

const vec2d vec2d::turned (const double& Winkel) const
{
	vec2d d=vec2d(cos(Winkel), -sin(Winkel));
	return (vec2d(*this*d,*this*!d));
}

const vec2d vec2d::unit_vector () const
{
	if (sqr(*this)<=grenze) return vec2d(0.0,0.0);

	else return (*this/vabs(*this));
}

const double vec2d::Richtung () const
{
	vec2d u=unit_vector();
	switch ((int)sign(u.X()))
	{
	   case -1:	return (sign(u.Y())*M_PI-asin(u.Y()));
	   case 1:	return (asin(u.Y()));
	   case 0:	return (M_PI_2*sign(u.Y()));
	   default: exit(0);
	}
}


const vec2d& vec2d::operator+=(const vec2d& V)
{
	x+= V.x;
	y+= V.y;
	return *this;
}

const vec2d& vec2d::operator-=(const vec2d& V)
{
	x-= V.x;
	y-= V.y;
	return *this;
}

const vec2d& vec2d::operator*=(const double& a)
{
	x*= a;
	y*= a;
	return *this;
}

const vec2d& vec2d::operator/=(const double& a)
{
	x/= a;
	y/= a;
	return *this;
}
