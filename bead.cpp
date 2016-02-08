#include "bead.h"

bead::bead(vektor<double> x)
{
	ort=x;
	xtimes=0;
	ytimes=0;
	ztimes=0;


}

bead::bead(double x, double y, double z,int xtimes, int ytimes, int ztimes)
{
    ort=vektor<double>(x,y,z);
    this->xtimes=xtimes;
    this->ytimes=ytimes;
    this->ztimes=ztimes;
}

bead::bead()
{
	double x=XMAX/2;
	double y=YMAX/2;
	double z=ZMAX/2;
	ort=vektor<double>(x,y,z);


}

double bead::Abstand(bead Nachbar)
{
	return ort.Abstand(Nachbar.ort);

}

double bead::Kopplungsenergie(bead Nachbar)
{
	return K*pow((ort-Nachbar.ort).Betrag()-x0,2);
}

void bead::verschiebe(vektor<double> v)
{
	if(ort.x+v.x>XMAX) xtimes++;
	if(ort.x+v.x<0) xtimes--;
	if(ort.y+v.y>YMAX) ytimes++;
	if(ort.y+v.y<0) ytimes--;
	if(ort.z+v.z>ZMAX) ztimes++;
	if(ort.z+v.z<0) ztimes--;
	ort.x=fmod(ort.x+v.x+XMAX,XMAX);
	ort.y=fmod(ort.y+v.y+YMAX,YMAX);
	ort.z=fmod(ort.z+v.z+ZMAX,ZMAX);

}
