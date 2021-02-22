/***************************************************************************
vec3d.h-description
 -------------------
begin: Wed May 22 2002
copyright: (C) 2002 by Malte Schmick
email: schmick@mpi-dortmund.mpg.de
 ***************************************************************************/

/***************************************************************************
 * *
 * This program is free software; you can redistribute it and/or modify*
 * it under the terms of the GNU General Public License as published by*
 * the Free Software Foundation; either version 2 of the License, or *
 * (at your option) any later version. *
 * *
 ***************************************************************************/

#ifndef VEC3D_H
#define VEC3D_H

# include <iostream>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <limits>
# include "myTempl.h"
using namespace std;

/**3D Vektoren mit double
*@author Malte Schmick
*/

class vec3d {
	protected:
	double x, y, z;

public:
	vec3d( double =0.0, double =0.0, double =0.0);			// Default-Konstruktor
	vec3d( vec3d *);
	vec3d( const vec3d&);
	~vec3d( void);										// Destruktor

	double X() const;
	double Y() const;
	double Z() const;
	
	void print() const
	{
	 	printf("( %f : %f : %f )\n", x, y, z);
	}
	void printStr(ostream& strm = cout) const;
	void InStr(istream& strm = cin);

	friend vec3d operator *(const double&, const vec3d&); 	// Skalar*Vektor=Vektor
		   vec3d operator *(const double&) const;			// Vektor*Skalar=Vektor
		   vec3d operator /(const vec3d&) const;			// VektorXVektor=Vektor
		   vec3d operator /(const double&) const;			// Vektor/Skalar=Vektor
		   vec3d operator +(const vec3d&) const;			// Vektor+Vektor=Vektor
		   vec3d operator -(const vec3d&) const;			// Vektor-Vektor=Vektor
		   vec3d operator -();												// Negation
	const  vec3d unit_vector () const;						// Einheitsvektor
	void normalize();

	double operator *(const vec3d&) const;			// Vektor*Vektor=Skalar
	const  double SQR (void) const;							// sqr(Vektor) als Objekt.Mitglied
	friend double vabs(const vec3d&); // Betragsabfrage als globale Funktion
	friend double sqr (const vec3d&);// Quadrieren als globale Funktion

	const int operator ==(const vec3d&) const;
	const int operator !=(const vec3d&) const;
	const vec3d& operator +=(const vec3d&);				// Vektor+=Vektor
	const vec3d& operator -=(const vec3d&);				// Vektor-=Vektor
	const vec3d& operator *=(const double&);			// Vektor/=double
	const vec3d& operator /=(const double&);			// Vektor/=double
};
ostream& operator<<(ostream& strm, const vec3d& V);
istream& operator>>(istream& strm, vec3d& V);

#endif
