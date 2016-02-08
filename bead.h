#ifndef BEAD_H
#define BEAD_H

//#include "vektor.cpp"

using namespace std;

class bead
{
    public:
	int xtimes;
	int ytimes;
	int ztimes;
        vektor<double> ort;
        bead(vektor<double> x);
        bead(double x, double y, double z,int xtimes, int ytimes, int ztimes);
        bead();
	void verschiebe(vektor<double> v);
	double Abstand(bead Nachbar);
	double Kopplungsenergie(bead Nachbar);
	double K;
	double x0;
    protected:
    private:
};

#endif // BEAD_H
