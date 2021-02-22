#include <limits.h>

#ifndef MYTEMPL_H
#define MYTEMPL_H
#define MY_DBL_DIG 12

using namespace std;

const double grenze = pow(10.,-2*MY_DBL_DIG+2);

const double rand_norm = 1./((double)RAND_MAX+1.);

template <class T>
	T sign (const T& a)	
	{
		return (sqr(a)<=grenze ? (T)0 : (a < 0 ? (T)(-1) : (T)(1)));
	}
template <class T>
	T MIN (const T& a, const T& b)	
	{
		return (a<=b ? a : b);
	}
template <class T>
	T MAX (const T& a, const T& b)	
	{
		return (a>=b ? a : b);
	}

template <class T>
	T CBC (const T& a)	
	{
		return (a*a*a);
	}


template <class T>
	T SQUARE (const T& a)	
	{
		return (a*a);
	}


#endif
