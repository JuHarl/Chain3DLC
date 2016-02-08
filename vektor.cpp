#include "vektor.h"

template <typename T> vektor<T>::vektor(T x, T y, T z)
{
	this->x=x;
	this->y=y;
	this->z=z;
}

template <typename T> void vektor<T>::Ausgabe()
{
	cout << "[" << x << "," << y << "," << z << "] \n" ;

}

template <typename T> T vektor<T>::Abstand(vektor<T> Abstandsvektor)//Berechnet den Abstand OHNE Berücksichtigung von periodischen Randbedingungen
{

	T x_=pow(x-Abstandsvektor.x,2);
	T y_=pow(y-Abstandsvektor.y,2);
	T z_=pow(z-Abstandsvektor.z,2);
	return (x_+y_+z_);

}

template <typename T> T vektor<T>::Betragsquadrat()
{
	return (x*x+y*y+z*z);

}

template <typename T> T vektor<T>::Betrag()
{
	return sqrt(x*x+y*y+z*z);

}

template <typename T> void vektor<T>::diff(vektor<T> addition)
{

	x-=addition.x;
	y-=addition.y;
	z-=addition.z;


    if(x< -XMAX/2) x+=XMAX;
    else if(x>XMAX/2) x-=XMAX;
    if(y< -YMAX/2) y+=YMAX;
    else if(y>XMAX/2) y-=YMAX;
    if(z< -ZMAX/2) z+=ZMAX;
    else if(z>ZMAX/2) z-=ZMAX;
//	if(abs(x_)>XMAX/2)
//	{
//		if(x_<0) x_+=XMAX;
//		else x_-=XMAX;
//	}
//	if(abs(y_)>YMAX/2)
//	{
//		if(y_<0) y_+=YMAX;
//		else y_-=YMAX;
//	}
//	if(abs(z_)>ZMAX/2)
//	{
//		if(z_<0) z_+=ZMAX;
//		else z_-=ZMAX;
//	}

//	vektor<T> neu=vektor<T>(x_, y_, z_);

//	return neu;
}
template <typename T> vektor<T> vektor<T>::operator-(vektor<T> addition)
{

	T x_=x-addition.x;
	T y_=y-addition.y;
	T z_=z-addition.z;


    if(x_< -XMAX/2) x_+=XMAX;
    else if(x_>XMAX/2) x_-=XMAX;
    if(y_< -YMAX/2) y_+=YMAX;
    else if(y_>XMAX/2) y_-=YMAX;
    if(z_< -ZMAX/2) z_+=ZMAX;
    else if(z_>ZMAX/2) z_-=ZMAX;
//	if(abs(x_)>XMAX/2)
//	{
//		if(x_<0) x_+=XMAX;
//		else x_-=XMAX;
//	}
//	if(abs(y_)>YMAX/2)
//	{
//		if(y_<0) y_+=YMAX;
//		else y_-=YMAX;
//	}
//	if(abs(z_)>ZMAX/2)
//	{
//		if(z_<0) z_+=ZMAX;
//		else z_-=ZMAX;
//	}

	vektor<T> neu=vektor<T>(x_, y_, z_);

	return neu;
}


template <typename T> vektor<T> vektor<T>::operator+(vektor<T> addition)
{

	T x_=x+addition.x;
	T y_=y+addition.y;
	T z_=z+addition.z;

	if(abs(x_)>XMAX)
	{
		if(x_<0) x_+=XMAX;
		else x_-=XMAX;
	}
	if(abs(y_)>YMAX)
	{
		if(y_<0) y_+=YMAX;
		else y_-=YMAX;
	}
	if(abs(z_)>ZMAX)
	{
		if(z_<0) z_+=ZMAX;
		else z_-=ZMAX;
	}

	vektor<T> neu=vektor<T>(x_, y_, z_);

	return neu;
}

template <typename T> vektor<T> vektor<T>::add(vektor<T> addition)//add ohne Berücksichtigung der periodischen Randbedingungen
{

	T x_=x+addition.x;
	T y_=y+addition.y;
	T z_=z+addition.z;


	vektor<T> neu=vektor<T>(x_, y_, z_);

	return neu;
}

template <typename T> vektor<T> vektor<T>::sub(vektor<T> addition)//sub ohne Berücksichtigung der periodischen Randbedingungen
{

	T x_=x-addition.x;
	T y_=y-addition.y;
	T z_=z-addition.z;


	vektor<T> neu=vektor<T>(x_, y_, z_);

	return neu;
}

template <typename T> T vektor<T>::operator*(vektor<T> addition)
{


	return x*addition.x+y*addition.y+z*addition.z;
}

template <typename T> vektor<T> vektor<T>::operator/(T div)
{

	T x_=x/div;
	T y_=y/div;
	T z_=z/div;
	vektor<T> neu=vektor<T>(x_, y_, z_);

	return neu;;
}

template <typename T> vektor<T> vektor<T>::operator*(T div)
{

	T x_=x*div;
	T y_=y*div;
	T z_=z*div;
	vektor<T> neu=vektor<T>(x_, y_, z_);

	return neu;
}

template <typename T> vektor<T> vektor<T>::kreuz(vektor<T> kreuz)
{
	T x_=y*kreuz.z-z*kreuz.y;
	T y_=z*kreuz.x-x*kreuz.z;
	T z_=x*kreuz.y-y*kreuz.x;

	vektor<T> neu=vektor<T>(x_,y_,z_);

	return neu;

}

template <typename T> void vektor<T>::verschiebe(vektor<T> v)
{
	x=fmod(x+v.x+XMAX,XMAX);
	y=fmod(y+v.y+YMAX,YMAX);
	z=fmod(z+v.z+ZMAX,ZMAX);

}
