#ifndef VEKTOR_H
#define VEKTOR_H



using namespace std;

template <typename T> class vektor
{
    public:
	T x;
	T y;
	T z;
	vektor(T x=0, T y=0, T z=0);
	//vektor();
	vektor add(vektor add);
	vektor sub(vektor add);
	vektor kreuz(vektor kreuz);
	T Betrag();
	T Abstand(vektor Abstandsvektor);
	void Ausgabe();
	void verschiebe(vektor v);
	void diff(vektor);
	vektor operator+(vektor);
	vektor operator-(vektor);
	vektor operator/(T);
	T Betragsquadrat();
	T operator*(vektor);//Skalarprodukt
	vektor operator*(T);

    protected:
    private:
};

#endif // VEKTOR_H
