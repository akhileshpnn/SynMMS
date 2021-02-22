#ifndef VEC2D_H
#define VEC2D_H


# include <iostream>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <limits.h>
# include "myTempl.h"

using namespace std;

class vec2d {

protected:
	double x,y;

public:
	vec2d(double =0.0, double =0.0);
	vec2d(vec2d *);
	vec2d(const vec2d&);
	~vec2d();

	double X() const;
	double Y() const;

	void print() const
	{
	 	printf("( %f : %f )\n", x, y);
	}
	void printStr(ostream& strm = cout) const;
	void InStr(istream& strm = cin);

	friend  vec2d  operator *(const double&, const vec2d&); 	// Skalar*Vektor=Vektor
	vec2d  operator *(const double&) const;			// Vektor*Skalar=Vektor
	vec2d  operator /(const double&) const;			// Vektor/Skalar=Vektor
	vec2d  operator +(const vec2d&) const;			// Vektor+Vektor=Vektor
	vec2d  operator -(const vec2d&) const;			// Vektor-Vektor=Vektor
	vec2d  operator -();												// Negation
	vec2d  operator !();												// Linksdrehung um 90grad
	const vec2d turned (const double&) const;		// Linksdrehung mit bel. Winkel
	const vec2d unit_vector () const;						// Einheitsvektor
	const double Richtung () const;							// Winkel zur x-Achse

	double operator *(const vec2d&) const;			// Vektor*Vektor=Skalar
	const  double SQR (void) const;							// sqr(Vektor) als Objekt.Mitglied
	friend double vabs (const vec2d&);         // Betragsabfrage als globale Funktion
	friend double sqr (const vec2d&);          // Quadrieren als globale Funktion

	const bool operator ==(const vec2d&) const;
	const bool operator !=(const vec2d&) const;
	const bool operator  >(const vec2d&) const;
	const bool operator  <(const vec2d&) const;
	const bool operator >=(const vec2d&) const;
	const bool operator <=(const vec2d&) const;
	const vec2d& operator +=(const vec2d&);			// Vektor+=Vektor
	const vec2d& operator -=(const vec2d&);			// Vektor-=Vektor
	const vec2d& operator *=(const double&);		// Vektor/=double
	const vec2d& operator /=(const double&);		// Vektor/=double
};

//Globale Überladung der Stream-Pipe-Operatoren
ostream& operator<<(ostream& strm, const vec2d& V);
istream& operator>>(istream& strm, vec2d& V);
#endif
