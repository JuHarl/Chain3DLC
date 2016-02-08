#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <functional>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <random>
#include <time.h>
#include <chrono>
#include <stdio.h>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "runningstat.cpp"
#include <thread>


#define _USE_MATH_DEFINES_

using namespace std;
using namespace Eigen;

double XMAX;//=10;
double YMAX;//=10;
double ZMAX;//=10;

const double Beta=1;
const double x0=1;
const double K=100;
long long tMax=1e6;
int Schrittweite=1e3;
int Sweeps=1e3;

const double l=2;

double rho=1;
int N=100;
int M;//rho*XMAX*YMAX*ZMAX/N;


int AnzZellen;//.8*XMAX;//1.0*XMAX/(1.25*x0);
double Zellgroesse;//1.0*XMAX/AnzZellen;


long long AnzahlKollisionen;
long long AnzahlFederEC;
long long AnzahlEndlosschleifen;

long long Schritt=0;

long long Bug=0;



const double PI=3.14159265359; //evtl noch aus Paket ziehen

mt19937 engine;
uniform_real_distribution<> distribution(0, 1);
auto rnd = bind(distribution, ref(engine));

uniform_real_distribution<> PIdistribution(0, PI);
auto rndPI = bind(PIdistribution, ref(engine));

uniform_real_distribution<> NEGdistribution(-1, 1);
auto rndWAY = bind(NEGdistribution, ref(engine));


//MTRand * Random;

#include "vektor.cpp"

int ** initMatrix(int Groesse)
{
    int** A = new int*[M];
     for(int i = 0; i < M; ++i)
        A[i] = new int[N];
    for(int i=0;i<M;i++){
        for(int j=0;j<N;j++){
            A[i][j]=0;
        }
    }
    return A;

}

double Federenergie(vektor<double> x)
{
	return 0.5*K*pow(x.Betrag()-x0,2);

}


double LGS(Matrix3f E, Vector3f G)
{
//	if(Schritt==Bug) cout << "Here is the matrix E:\n" << E << endl;
//	if(Schritt==Bug) cout << "Here is the vector G:\n" << G << endl;

	FullPivHouseholderQR<Matrix3f> dec(E);
	Vector3f x = dec.solve(G);

	//if(Schritt==Bug) cout << "Lösungsvektor x=\n "<< x << endl;

  	double relative_error = (E*x - G).norm() / G.norm(); // norm() is L2 norm
//   	if(Schritt==Bug) cout << "The relative error is:\n" << relative_error << endl;

	if(relative_error>1e-4) return 1000;
	if(abs(x[2])/relative_error<1) return 0;
	return x[2];
}

double SchnittpunktEbeneGerade(vektor<double> b0, vektor<double> b1, vektor<double> r, vektor<double> ev)//f ist Nachbar-bead von r, b0,b1 sind Verbindungsvektoren von r zu dem kollisionsbeads
{
	Matrix3f E;
	E << b0.x,b1.x,-ev.x,	b0.y,b1.y,-ev.y,	b0.z,b1.z,-ev.z;
	Vector3f G;
	G << r.x,r.y,r.z;



	FullPivHouseholderQR<Matrix3f> dec(E);
	Vector3f x = dec.solve(G);

	//if(Schritt==Bug) cout << "Lösungsvektor x=\n "<< x << endl;

  	double relative_error = (E*x - G).norm() / G.norm(); // norm() is L2 norm
//   	if(Schritt==Bug) cout << "The relative error is:\n" << relative_error << endl;

	if(relative_error>1e-4) return 1000;
	if(abs(x[2])/relative_error<1) return 0;
	return x[2];
}

double LGS(vektor<double> b0, vektor<double> b1, vektor<double> r, vektor<double> ev)
	{
		double detA=(b0.z*b1.y*ev.x-b0.y*b1.z*ev.x-b0.z*b1.x*ev.y+b0.x*b1.z*ev.y+b0.y*b1.x*ev.z-b0.x*b1.y*ev.z);
		if(detA==0) return 1000;
		return -(b0.z*b1.y*r.x-b0.y*b1.z*r.x-b0.z*b1.x*r.y+b0.x*b1.z*r.y+b0.y*b1.x*r.z-b0.x*b1.y*r.z)/detA; //Cramersche Regel
	}

pair<double,double> ProjektionEbene(vektor<double> b0, vektor<double> b1, vektor<double> r) //Projiziert r auf Ebene aufgespannt durch b0, b1
{
	Matrix2f E;
	Vector2f G;
	E << b0*b0, b1*b0,	b0*b1, b1*b1;
	G << r*b0, r*b1;

	FullPivHouseholderQR<Matrix2f> dec(E);
	Vector2f x = dec.solve(G);

  	//double relative_error = (E*x - G).norm() / G.norm(); // norm() is L2 norm
   	//if(Schritt==Bug) cout << "The relative error is:\n" << relative_error << endl;

	//if(relative_error>1e-4) return 1000;
	//if(abs(x[2])/relative_error<10) return 0;
	return make_pair(x[0],x[1]);

}


/*
pair<double,double> SchnittpunktGeradeGerade(vektor<double> S, vektor<double> b0, vektor <double> b1)
{
	S.Ausgabe();
	b0.Ausgabe();
	b1.Ausgabe();
	double q,p;
	if(b1.y!=b0.y){ //ACHTUNG!
		if(S.x==0){
			p=b0.y/(b0.y-b1.y);
			if(S.y!=0)
			{
				q=(b0.y+(b1.y-b0.y)*p)/S.y;
			}
			else {
				q=(b0.z+(b1.z-b0.z)*p)/S.z;
			}

		}
		else{
			double B=(S.y)/(S.x);
			double Q=(b0.y)/(b1.y-b0.y)-B*(b0.x)/(b1.y-b0.y);
			p=-Q/(1-B*(b1.x-b0.x)/(b1.y-b0.y));
			q=(b0.x+(b1.x-b0.x)*p)/(S.x);
			}
		}
	else
	{	if(S.x==0){
			p=b0.x/(b0.x-b1.x);
			if(S.y!=0)
			{
				q=(b0.y+(b1.y-b0.y)*p)/S.y;
			}
			else {
				q=(b0.z+(b1.z-b0.z)*p)/S.z;
			}

		}
			else{
		double B=(S.z)/(S.x);
		double Q=(b0.z)/(b1.z-b0.z)-B*(b0.x)/(b1.z-b0.z);
		p=-Q/(1-B*(b1.x-b0.x)/(b1.z-b0.z));
		q=(b0.x+(b1.x-b0.x)*p)/(S.x);
		}
	}
	//if(abs(p)==nan || abs(q)==nan) cout << "FEHLER!!!!!!!" << endl;
	cout << "p=" << p << "\t q="<< q << endl;
	return make_pair(p,q);
}*/

double SchnittpunktGeradeGerade(vektor<double> S, vektor<double> b0, vektor <double> b1)
{
//	S.Ausgabe();
//	b0.Ausgabe();
//	b1.Ausgabe();
	double Sb0=(S*b0)/(S.Betrag()*b0.Betrag());
	double Sb1=(S*b1)/(S.Betrag()*b1.Betrag());
	double b0b1=(b0*b1)/(b0.Betrag()*b1.Betrag());
	double phi0=acos(Sb0);
	double phi1=acos(Sb1);
	double phi2=acos(b0b1);
//	if(Schritt==Bug) cout << "phi0=" << phi0/PI*180 << "° phi1=" << phi1/PI*180 << "° phi2=" << phi2/PI*180 << "°" << endl;
//	if(Schritt==Bug) cout << "Sb0=" << Sb0 << " Sb1=" << Sb1 << " b0b1=" << b0b1 << endl;
	if(abs(phi2)<1.e-7 || abs(b0b1-1.)<1.e-7 ) return 1000.;
	if(abs(Sb1-1.)<1.e-7) return 1.;
	if(abs(Sb0-1.)<1.e-7) return 0.;
	if(Sb0<0 && Sb1<0) return -1; //Schnittpunkt liegt entgegen des Bonds
	if(abs(phi0+phi1+phi2-2*PI)<1e-5) return -1.5;
	if(phi1/phi2>1.) return 2.;
	if(phi0/phi2>1.) return 2.;

	double gamma=acos((-(b0*(b1-b0)))/(b0.Betrag()*(b1-b0).Betrag()));
	double delta=PI-phi0-gamma;
	double p=b0.Betrag()/sin(delta)*sin(phi0);
	return p/(b1-b0).Betrag();

}


pair<double,double> Kollision(vektor<double> f, vektor<double> b0, vektor<double> b1, vektor<double> r, vektor<double> ev,double x) //1. Parameter: distance, 2. Parameter: p
{
//	if(Schritt==Bug){
//		cout << "Vor Subtraktion von f: " << endl;
//		cout << "b0=";
//		b0.Ausgabe();
//		cout << "b1=";
//		b1.Ausgabe();
//		cout << "r=";
//		r.Ausgabe();
//		cout << "f=";
//		f.Ausgabe();
//		cout << "ev=";
//		ev.Ausgabe();
//	}
	vektor<double> b=b0-b1;

	b1.diff(f);
	b0=b1.add(b);
//	r=r-f;
//	bool ausgeschlossen=false;

	if(ev.x>0){
		if(b0.x<r.x && b0.x<0 && b1.x<r.x && b1.x<0) {
//			if(Schritt==Bug) cout << "x1" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
		}
		if(b0.x>r.x+ev.x*x && b0.x>ev.x*x && b1.x>r.x+ev.x*x && b1.x>ev.x*x)
		{
//			if(Schritt==Bug) cout << "x2" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
		}
	}
	else{
		if(b0.x<r.x+ev.x*x && b0.x < ev.x*x && b1.x<r.x+ev.x*x && b1.x < ev.x) {
//			if(Schritt==Bug) cout << "x3" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
			}
		if(b0.x>r.x && b0.x>0 && b1.x>r.x&&b1.x>0) {
//			if(Schritt==Bug) cout << "x4" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
			}

	}
	if(ev.y>0){
		if(b0.y<r.y && b0.y<0 && b1.y<r.y && b1.y<0){
//			if(Schritt==Bug) cout << "y1" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
			}
		if(b0.y>r.y+ev.y*x && b0.y>ev.y*x && b1.y>r.y+ev.y*x && b1.y>ev.y*x){
//			if(Schritt==Bug)cout << "y2" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
			}
	}
	else{
		if(b0.y<r.y+ev.y*x && b0.y < ev.y*x && b1.y<r.y+ev.y*x && b1.y < ev.y)
			{
//			if(Schritt==Bug) cout << "y3" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
			}
		if(b0.y>r.y && b0.y>0 && b1.y>r.y&&b1.y>0) {
//			if(Schritt==Bug) cout << "y4" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
			}
	}
	if(ev.z>0){
		if(b0.z<r.z && b0.z <0 && b1.z<r.z && b1.z<0) {
//			if(Schritt==Bug) cout << "z1" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
			}
		if(b0.z>r.z+ev.z*x && b0.z>ev.z*x && b1.z>r.z+ev.z*x && b1.z>ev.z*x) {
//			if(Schritt==Bug) cout << "z2" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
			}
	}
	else{
		if(b0.z<r.z+ev.z*x && b0.z < ev.z*x && b1.z<r.z+ev.z*x && b1.z < ev.z) {
//			if(Schritt==Bug) cout << "z3" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
			}
		if(b0.z>r.z && b0.z>0 && b1.z>r.z&&b1.z>0){
//			if(Schritt==Bug) cout << "z4" << endl;
//			ausgeschlossen=true;
			return make_pair(x,10);
			}

	}

//	if(Schritt==Bug){
//		cout << "Nach Subtraktion von f: " << endl;
//		cout << "b0=";
//		b0.Ausgabe();
//		cout << "b1=";
//		b1.Ausgabe();
//		cout << "r=";
//		r.Ausgabe();
//	}
	double w=LGS(b0,b1,r,ev);
//	if(Schritt==Bug) cout << "w= " << w << endl;
//	double w=SchnittpunktEbeneGerade(b0,b1,r,ev);
//	if(Schritt==Bug) cout << "w= " << w << endl;
	if(w>=x)
		{
//			if(Schritt==Bug) cout << " Schnittpunkt so weit weg, dass keine Kollision möglich " << endl;
			return make_pair(x,10); //Schnittpunkt so weit weg, dass keine Kollision möglich
		}
	if(w<=0)//Toleranz falls direkt aneinander --- MUSS NOCH ÜBERPRÜFT WERDEN!!!
		{

//			if(Schritt==Bug) cout << "Bewegung in falsche Richtung, keine Kollision möglich " << endl;
			return make_pair(x,10);
		}
	double p,q;
	vektor<double> SG;
//	if(w==0) {//r in Ebene
////		if(Schritt==Bug) cout << "ACHTUNG! w ist exakt 0" << endl;
//		if(abs((b0.kreuz(b1))*ev)==0){//ev ebenfalls in Ebene
////			if(Schritt==Bug) {
////				b0.kreuz(b1).Ausgabe();
////				cout << "(b0 x b1) * ev =" << b0.kreuz(b1)*ev << endl;
////				cout << "ev ebenfalls in Ebene " << endl;
////			}
//			/*if(ev*b0.sub(r)<0 && ev*b1.sub(r)<0) {
//				cout << "Bewegung von bond weg" << endl;
//				return x;
//
//			}*/
//			p=SchnittpunktGeradeGerade(r.add(ev*x),b0,b1);
////			if(Schritt==Bug) cout << "p=" << p << endl;
//			if(!(p<=1) || !(p>=0)) {//AUFPASSEN!!! Fängt p=nan mit ab!
////				if(Schritt==Bug) cout << "Bead bewegt sich an Kollisionsbond vorbei" << endl;
//				return x; //Bead bewegt sich an Kollisionsbond vorbei
//				}
//
//			SG=b0*(1.-p)+b1*p;
////			if(Schritt==Bug) SG.Ausgabe();
//			w=SG.sub(r).Betrag();
////			if(Schritt==Bug) cout << "w=" << w << endl;
//			if(w<0) return x;
////			if(Schritt==Bug && w<x) cout << "\t !!!!!!!!!!!!!!!!!!!!!!!!!Kollision in Ebene!!!!!!!!!!!!!!!!" << endl;
//			return min(w,x);
//		}
//		else{
////			if(Schritt==Bug) {
////				//b0.kreuz(b1).Ausgabe();
////				//cout << "(b0 x b1) * ev =" << b0.kreuz(b1)*ev << endl;
////				cout << "ev nicht in Ebene , keine Kollision möglich" << endl;
////			}
//			return x;
//			}
//		}
	vektor<double> S=r.add(ev*w);
//	if(Schritt==Bug) S.Ausgabe();
	//double k0,k1;
	//pair<double,double> k0k1=ProjektionEbene(b0,b1,S);
	//k0=k0k1.first;
	//k1=k0k1.second;
	//S=b0*k0+b1*k1;
	//if(Schritt==Bug) S.Ausgabe();
	//pair<double,double> pq=SchnittpunktGeradeGerade(S,b0,b1);
	//p=pq.first;
	//q=pq.second;
	p=SchnittpunktGeradeGerade(S,b0,b1);
//	if(Schritt==Bug) cout << "p=" << p << endl;

	//double p_=(int)(p*1000+.5)/1000.0; //Rundet p auf 3 Nachkommastellen... ACHTUNG!
	//if(p==1) cout << "p==1" << endl;

	if(!(p<=1.+1.e-5) || !(p>=0.-1.e-5)) {//AUFPASSEN!!! Fängt p=nan mit ab!
//			if(Schritt==Bug) cout << "Bead bewegt sich an Kollisionsbond vorbei" << endl;
			return make_pair(x,10); //Bead bewegt sich an Kollisionsbond vorbei
			}
	/*if(q<0)	{
			cout << "Schnittpunkt hinter f, keine Kollision" << endl;
			return x;
		}
	*/
	//p=.5;
	SG=b0*(1.-p)+b1*p;
// 	if(Schritt==Bug) SG.Ausgabe();
//	if(Schritt==Bug) cout << "(SG).Betrag()=" << (SG).Betrag() << " (S).Betrag()=" << (S).Betrag() << endl;
	if((SG).Betrag()>(S).Betrag()+1.e-5) {//Toleranz von 1e-2!! ACHTUNG!

//		if(Schritt==Bug) cout << "Bead bewegt sich vor dem Bond vorbei" << endl;
		return make_pair(x,10);//Bead bewegt sich vor dem Bond vorbei
		}
	else {
//	if(Schritt==Bug) cout << "\t !!!!!!!!!!!!!!!!!!!!!!!!!Kollision!!!!!!!!!!!!!!!!" << endl;
//	if(ausgeschlossen) cout << "FEHLER! KOLLISION MIT AUSGESCHLOSSENEM BOND!" << endl;
	return make_pair(w,p);//Kollision
	}
}

//class chain;
//int k;

#include "bead.cpp"

#include "chain.h"

//vector<chain> Chains;


void Positionen(ofstream &ofbla, vector<chain> * Chains);
#include "chain.cpp"

#include "contactmap.cpp"


void Positionen(ofstream &ofbla, vector<chain> * Chains){

	ofbla <<"-------------A-N-F-A-N-G--------------------" << endl << endl;

        ofbla << "Datum: " << 1 << "." << 11 << "." << 11 <<"\t" << "Zeit:  "  << 1 <<
":" << 1 << ":" << 1 << endl;
        ofbla << "Anz: " << M << "\tN: " << N-1 << "\tb0: " << x0  <<
"\tFeder: " << K << "\tBiege: " << 0 << "\tPotential: "
<< 0 << "\tTiefe: " << 0 << "\thard: " << 1 <<
"\tBeta: " << Beta << "\tWand: " << 0 << "\tL_1: " <<
XMAX << "\tL_2: " << YMAX << "\tL_3(leer): " << ZMAX  << "\tleer: "
<< "0"<< "\tleer: " << "0"<< "\tleer: " << "0"  << endl<<endl;
	for(int j=0;j<M;j++){
		    for(int i =0;i<N;i++)
		    {
		       ofbla << Chains->at(j).Beads.at(i).ort.x << "\t" << Chains->at(j).Beads.at(i).ort.y <<
	"\t"<< Chains->at(j).Beads.at(i).ort.z << endl;
	//               cout << k << " " << i << " " << j << endl;
		    }
		    ofbla << endl;
		}

        ofbla <<"---------------E-N-D-E----------------------" << endl << endl;

}

void g1(ofstream &o, vector<vektor<double>> MidBeads)
{
    cout.precision(100);
    cout << "Berechne g1...";
    RunningStat RSg1;
    for(int i=1;i<Schritt;i*=sqrt(sqrt(2)))
    {
        for(int j=0;j<M;j++)
        {
            for(int k=i;k<Schritt;k+=i){
                RSg1.Push(MidBeads.at(k*M+j).Abstand(MidBeads.at((k-i)*M+j)));
//                    cout << i << "\t" << j << "\t" << k << endl;
//                if(i==1) cout << MidBeads.at(k*M+j).Abstand(MidBeads.at((k-i)*M+j)) << endl;
//                    if(i==1 && j==M-1 && k==Schritt-1){
//                    cout << MidBeads.at(k*M+j).x << "\t" << MidBeads.at(k*M+j).y << "\t" << MidBeads.at(k*M+j).z << "\t" << endl;
//                    cout << MidBeads.at((k-i)*M+j).x << "\t" << MidBeads.at((k-i)*M+j).y << "\t" << MidBeads.at((k-i)*M+j).z << "\t" << endl;
//
//                }
            }
        }
        o << i << "\t" << RSg1.Mean() << "\t"<< RSg1.StandardDeviation() << "\n";
        //cout << "g1 gemittelt über " << RSg1.NumDataValues() << "für t=" << i << endl;
        RSg1.Clear();
        i++;
    }
    cout << "fertig" << endl;


}
void g2(ofstream &o, vector<vektor<double>> MidBeads, vector<vektor<double>> Schwerpunkte)
{
    cout << "Berechne g2...";


    RunningStat RSg1;
    for(int i=1;i<Schritt;i*=sqrt(sqrt(2)))
    {
        for(int j=0;j<M;j++)
        {
            for(int k=i;k<Schritt;k+=i){
                RSg1.Push(MidBeads.at(k*M+j).sub(Schwerpunkte.at(k*M+j)).Abstand(MidBeads.at((k-i)*M+j).sub(Schwerpunkte.at((k-i)*M+j))));
            }
        }
        o << i << "\t" << RSg1.Mean() << "\t"<< RSg1.StandardDeviation() << "\n";
//        cout << "g1 gemittelt über " << RSg1.NumDataValues() << "für t=" << i << endl;
        RSg1.Clear();
        i++;
    }
    cout << "fertig" << endl;



}

void g3(ofstream &o, vector<vektor<double>> Schwerpunkte)
{
    cout << "Berechne g3...";

    RunningStat RSg1;
    for(int i=1;i<Schritt;i*=sqrt(sqrt(2)))
    {
        for(int j=0;j<M;j++)
        {
            for(int k=i;k<Schritt;k+=i){
                RSg1.Push(Schwerpunkte.at((k)*M+j).Abstand(Schwerpunkte.at((k-i)*M+j)));
            }
        }
        o << i << "\t" << RSg1.Mean() << "\t"<< RSg1.StandardDeviation() << "\n";
//        cout << "g1 gemittelt über " << RSg1.NumDataValues() << "für t=" << i << endl;
        RSg1.Clear();
        i++;
    }
    cout << "fertig" << endl;



}
void g4(ofstream &o, vector<vektor<double>> FirstBeads, vector<vektor<double>> LastBeads)
{
    cout << "Berechne g4...";

    RunningStat RSg1;
    for(int i=1;i<Schritt;i*=sqrt(sqrt(2)))
    {
        for(int j=0;j<M;j++)
        {
            for(int k=i;k<Schritt;k+=i){
                RSg1.Push(FirstBeads.at((k)*M+j).Abstand(FirstBeads.at((k-i)*M+j)));
                RSg1.Push(LastBeads.at((k)*M+j).Abstand(LastBeads.at((k-i)*M+j)));
            }
        }
        o << i << "\t" << RSg1.Mean() << "\t"<< RSg1.StandardDeviation() << "\n";
//        cout << "g1 gemittelt über " << RSg1.NumDataValues() << "für t=" << i << endl;
        RSg1.Clear();
        i++;
    }
    cout << "fertig" << endl;



}

void g5(ofstream &o, vector<vektor<double>> FirstBeads, vector<vektor<double>> LastBeads, vector<vektor<double>> Schwerpunkte)
{
    cout << "Berechne g5...";

    RunningStat RSg1;
    for(int i=1;i<Schritt;i*=sqrt(sqrt(2)))
    {
        for(int j=0;j<M;j++)
        {
            for(int k=i;k<Schritt;k+=i){
                RSg1.Push(FirstBeads.at((k)*M+j).sub(Schwerpunkte.at((k)*M+j)).Abstand(FirstBeads.at((k-i)*M+j).sub(Schwerpunkte.at((k-i)*M+j))));
                RSg1.Push(LastBeads.at((k)*M+j).sub(Schwerpunkte.at((k)*M+j)).Abstand(LastBeads.at((k-i)*M+j).sub(Schwerpunkte.at((k-i)*M+j))));
            }
        }
        o << i << "\t" << RSg1.Mean() << "\t"<< RSg1.StandardDeviation() << "\n";
//        cout << "g1 gemittelt über " << RSg1.NumDataValues() << "für t=" << i << endl;
        RSg1.Clear();
        i++;
    }
    cout << "fertig" << endl;



}

void g6(ofstream &o, vector<vektor<double>> MidBeads)
{
    cout << "Berechne g6...";

    RunningStat RSg1;
    for(int i=1;i<Schritt;i*=sqrt(sqrt(2)))
    {
        for(int j=0;j<M;j++)
        {
            for(int k=i;k<Schritt;k+=i){
                    double x=0;
                    x+=pow(MidBeads.at((k)*M+j).x-MidBeads.at((k-i)*M+j).x,4);
                    x+=pow(MidBeads.at((k)*M+j).y-MidBeads.at((k-i)*M+j).y,4);
                    x+=pow(MidBeads.at((k)*M+j).z-MidBeads.at((k-i)*M+j).z,4);
                    RSg1.Push(x);
            }
        }
        o << i << "\t" << sqrt(RSg1.Mean()) << "\t"<< sqrt(RSg1.StandardDeviation()) << "\n";
//        cout << "g1 gemittelt über " << RSg1.NumDataValues() << "für t=" << i << endl;
        RSg1.Clear();
        i++;
    }
    cout << "fertig" << endl;



}

void BerechneGs(vector<chain> *Chains)
{
    cout.precision(50);
    ofstream o;
    ifstream datei("MidBeads.txt",ios::in);


    vector<vektor<double>> Schwerpunkte;
    vector<vektor<double>> MidBeads;
    vector<vektor<double>> FirstBeads;
    vector<vektor<double>> LastBeads;
    double x;
    double y;
    double z;

    while(!datei.eof()){
        datei >> x;
        datei >> y;
        datei >> z;
        Schwerpunkte.push_back(vektor<double>(x,y,z));
        datei >> x;
        datei >> y;
        datei >> z;
//        cout << "MidBead.push_back:" << x << "\t" << y << "\t" << z << "\n";
        MidBeads.push_back(vektor<double>(x,y,z));
        datei >> x;
        datei >> y;
        datei >> z;
        FirstBeads.push_back(vektor<double>(x,y,z));
        datei >> x;
        datei >> y;
        datei >> z;
        LastBeads.push_back(vektor<double>(x,y,z));

    }
    o.open("g1.txt");
    g1(o,MidBeads);
    o.close();
    o.open("g2.txt");
    g2(o,MidBeads,Schwerpunkte);
    o.close();
    o.open("g3.txt");
    g3(o,Schwerpunkte);
    o.close();
    o.open("g4.txt");
    g4(o,FirstBeads,LastBeads);
    o.close();
    o.open("g5.txt");
    g5(o,FirstBeads,LastBeads,Schwerpunkte);
    o.close();
    o.open("g6.txt");
    g6(o,MidBeads);
    o.close();
}
/*

void MidBeadKorr(ofstream &o, vector<chain> * Chains)
{
	for(int i=0;i<M;i++) Chains->at(i).MidBeadKorr();//Midbeadkorrelationen der einzelnen Chains ausrechnen
	RunningStat RS;
	int t=1;
	for(int i=0;i<Chains->at(0).midBeadKorrMeans.size() ;i++){
		for(int j=0;j<M;j++)
			{
			RS.Push(Chains->at(j).midBeadKorrMeans.at(i));
			}
		o << t << "\t" << RS.Mean() << "\t" << RS.StandardDeviation() << "\n" ;
//        cout << "MidBearKorr gemittelt über " << RS.NumDataValues() << "für t=" << t << endl;

		RS.Clear();
		t++;
		t*=sqrt(sqrt(2));
	}



}

*/
//vector<chain> Chains;

//void MessungEndtoEnd(int Nmin,int Nmax, int NSchritt)
//	{
//		RunningStat RSEndtoEnd;
//		ofstream oEE;
//		oEE.open("EndToEnd.txt");
//		vector<chain> * Chains = new vector<chain>();;
//		for(int k=Nmin;k<=Nmax;k+=NSchritt){
//			Chains->push_back(chain(k, Chains));
//			cout << "N=" << k << endl;
//			N=k;
//			for(int j=0;j<tMax;j+=Schrittweite)
//				{

//					double phi=2*rndPI();
//					double theta=rndPI();
//					vektor<double> ev=vektor<double>(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta));
//			//		vektor<double> ev=vektor<double>(.1,-1,0);
//					ev=ev/ev.Betrag();
//					int i=rnd()*k;
//					for(int q=0;q<Schrittweite;q++) Chains->at(0).EC(l,i,0,0,i,ev);
//					RSEndtoEnd.Push(Chains->at(0).EndToEnd());
//				}
//				oEE << k << "\t" << RSEndtoEnd.Mean() << "\t" << RSEndtoEnd.StandardDeviation() << "\n";
//				RSEndtoEnd.Clear();
//				Chains->clear();
//			}
//

//	}

void EndtoEndThread(int k, double * Array, int index)
	{

			mt19937 engine;

			uniform_real_distribution<> distribution(0, 1);
			auto rnd = bind(distribution, ref(engine));

			uniform_real_distribution<> PIdistribution(0, PI);
			auto rndPI = bind(PIdistribution, ref(engine));

			uniform_real_distribution<> NEGdistribution(-1, 1);
			auto rndWAY = bind(NEGdistribution, ref(engine));


//			cout << "Seed" << endl;
		   	time_t t;
		    	time(&t);
		    	srand((unsigned int)t);
		    	engine.seed(t+index*1000);


			vector<chain> * Chains = new vector<chain>();

            list<pair<int,int>> ***Zellen;

              // Allocate memory
              Zellen = new list<pair<int,int>>**[AnzZellen];
              for (int i = 0; i < AnzZellen; ++i) {
                Zellen[i] = new list<pair<int,int>>*[AnzZellen];

                for (int j = 0; j < AnzZellen; ++j)
                  Zellen[i][j] = new list<pair<int,int>>[AnzZellen];
              }

			Chains->push_back(chain(k,Chains,&engine,Zellen));


//			cout << "Startkonfiguration erzeugt " << endl;

//			ofstream ofbla;
//			ofbla.open("Pos.txt");
			for(int q=0;q<Schrittweite*k*tMax;q++){//Äquilibrierung
//						cout << "q=" << q << endl;
						Schritt=q;
						double phi=2*rndPI();
						double theta=acos(1-2*rnd());
						vektor<double> ev=vektor<double>(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta));
				//		vektor<double> ev=vektor<double>(.1,-1,0);
						ev=ev/ev.Betrag();
						int i=rnd()*k;
						Chains->at(0).EC(l,i,0,0,i,ev);
//						Positionen(ofbla,Chains);
				}

//			cout << "Äquilibrierung abgeschlossen" << endl;
			for(int j=0;j<tMax;j+=Schrittweite)
				{
					for(int q=0;q<Schrittweite*k;q++){
						double phi=2*rndPI();
						double theta=acos(1-2*rnd());
						vektor<double> ev=vektor<double>(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta));
				//		vektor<double> ev=vektor<double>(.1,-1,0);
						ev=ev/ev.Betrag();
						int i=rnd()*k;
						Chains->at(0).EC(l,i,0,0,i,ev);
					}
//					cout << "j=" << j << endl;
					Array[j/Schrittweite+index]=Chains->at(0).EndToEnd();
//					RSEndtoEnd->Push(Chains->at(0).EndToEnd());
//					if(Chains->at(0).EndToEnd()>k*k){
//					cout << "Schritt " << j << "EndToEnd=" << Chains->at(0).EndToEnd() << endl;
//					Chains->at(0).Beads.at(N-1).ort.Ausgabe();
//					Chains->at(0).Beads.at(0).ort.Ausgabe();

//					}

//					if(Chains->at(0).EndToEnd()>1000) {cout << "Schritt " << j << "EndToEnd=" << Chains->at(0).EndToEnd() << endl;

//					Chains->at(0).Beads.at(N-1).ort.Ausgabe();
//					Chains->at(0).Beads.at(0).ort.Ausgabe();	}
//					(Chains->at(0).Beads.at(N-1).ort-Chains->at(0).Beads.at(0).ort).Ausgabe();
//					cout << Chains->at(0).EndToEnd() << endl;
//					Chains->at(0)
//					if(Chains->at(0).EndToEnd()!=(Chains->at(0).Beads.at(N-1).ort-Chains->at(0).Beads.at(0).ort).Betragsquadrat()) cout << "Fehler bei Berechnung von EndToEnd" << endl;
				}
			delete Chains;
				  // De-Allocate memory to prevent memory leak
              for (int i = 0; i < AnzZellen; ++i) {
                for (int j = 0; j < AnzZellen; ++j)
                  delete [] Zellen[i][j];

                delete [] Zellen[i];
              }
              delete [] Zellen;
			cout << "Thread " << index/(tMax/Schrittweite)+1 << " beendet" << endl;

	}

void MessungEndtoEndThreads(int Nmin,int Nmax, int NSchritt,bool Verteilung=false)//M=1 setzen!
	{
		RunningStat * RSEndtoEnd =new RunningStat();
		ofstream oEE;
		ofstream oDist;
		oEE.open("EndToEnd.txt");
		double *EndToEnds = new double[tMax/Schrittweite*8];
		int AnzKerne=8;
		int AnzJ=1;


		for(int k=Nmin;k<=Nmax;k*=NSchritt){
				N=k;
				cout << "N=" << k << endl;
				if(Verteilung){
						oDist.open("VerteilungN_"+to_string(N)+".txt");
				}
				for(int j=0;j<AnzJ;j++){
					if(j%1==0) cout << "j=" << j << endl;
//					EndtoEndThread(k,EndToEnds,0);
					thread t1(EndtoEndThread,k,EndToEnds,0);
					thread t2(EndtoEndThread,k,EndToEnds,tMax/Schrittweite*1);
					thread t3(EndtoEndThread,k,EndToEnds,tMax/Schrittweite*2);
					thread t4(EndtoEndThread,k,EndToEnds,tMax/Schrittweite*3);
					thread t5(EndtoEndThread,k,EndToEnds,tMax/Schrittweite*4);
					thread t6(EndtoEndThread,k,EndToEnds,tMax/Schrittweite*5);
					thread t7(EndtoEndThread,k,EndToEnds,tMax/Schrittweite*6);
					thread t8(EndtoEndThread,k,EndToEnds,tMax/Schrittweite*7);

					t1.join();
					t2.join();
					t3.join();
					t4.join();
					t5.join();
					t6.join();
					t7.join();
					t8.join();

				for(int i=0;i<tMax/Schrittweite*AnzKerne;i++) {
					RSEndtoEnd->Push(EndToEnds[i]);
					if(Verteilung){
						oDist << EndToEnds[i] << endl;

					}
				}

				}
				if(Verteilung) oDist.close();
//				EndtoEndThread(k,RSEndtoEnd);
				cout << RSEndtoEnd->NumDataValues() << endl;

				oEE << k-1 << "\t" << RSEndtoEnd->Mean() << "\t" << RSEndtoEnd->StandardDeviation() << "\n";
				RSEndtoEnd->Clear();
//				Chains.clear();
//				cout << "Chains.clear()" << endl;
			}
		delete RSEndtoEnd;


	}

void LCOrientierung(ofstream &o, vector<chain> * Chains)
	{
		RunningStat RSO;
		vektor<double> O=vektor<double>(0,0,0);
		vektor<double> N=vektor<double>(0,0,0);
		for(int i=0;i<M;i++)
			{
				vektor<double> q=Chains->at(i).Beads.at(0).ort-Chains->at(i).Beads.at(1).ort;
				N=N.add(q);

			}
		N=N/N.Betrag();
		for(int i=0;i<M;i++)
			{
				vektor<double> q=Chains->at(i).Beads.at(0).ort-Chains->at(i).Beads.at(1).ort;
				q=q/q.Betrag();
				double Cosin=q*N;
				RSO.Push(Cosin*Cosin);
			}
		o << .5*(3*RSO.Mean()-1) << "\t" << 1.5*RSO.StandardDeviation() << "\n";
	}

void LCOrientierung2(ofstream &o, vector<chain> * Chains)
	{
		RunningStat RSO;
		Matrix3f Q = Matrix3f::Zero();
		for(int i=0;i<M;i++)
			{
				vektor<double> q=Chains->at(i).Beads.at(0).ort-Chains->at(i).Beads.at(1).ort;
				q=q/q.Betrag();
				Q(0,0)+=q.x*q.x-1./3;
				Q(0,1)+=q.x*q.y;
				Q(0,2)+=q.x*q.z;
				Q(1,1)+=q.y*q.y-1./3;
				Q(1,2)+=q.y*q.z;
				Q(2,2)+=q.z*q.z-1./3;

			}
		Q(1,0)=Q(0,1);
		Q(2,0)=Q(0,2);
		Q(2,1)=Q(1,2);
		Q/=M;
//		cout << "Q=" << Q << endl;
		SelfAdjointEigenSolver<Matrix3f> es(Q,false);
		double S=0;
		for(int i=0;i<3;i++){
				if(es.eigenvalues()(i)>S) S=es.eigenvalues()(i);
			}
		o << S*3./2 << "\t" << 0 << "\n" ;
	}

void MesseOrientierung(ofstream &o, vector<chain> * Chains, ofstream &ofbla,bool MessePos)
	{
		RunningStat RSO;
		Matrix3f Q = Matrix3f::Zero();
		for(int i=0;i<M;i++)
				{
					vektor<double> q=Chains->at(i).Beads.at(0).ort-Chains->at(i).Beads.at(1).ort;
					q=q/q.Betrag();
					Q(0,0)+=q.x*q.x-1./3;
					Q(0,1)+=q.x*q.y;
					Q(0,2)+=q.x*q.z;
					Q(1,1)+=q.y*q.y-1./3;
					Q(1,2)+=q.y*q.z;
					Q(2,2)+=q.z*q.z-1./3;

				}
			Q(1,0)=Q(0,1);
			Q(2,0)=Q(0,2);
			Q(2,1)=Q(1,2);
			Q/=M;
	//		cout << "Q=" << Q << endl;
			SelfAdjointEigenSolver<Matrix3f> es(Q,false);
			double S=0;
			for(int i=0;i<3;i++){
					if(es.eigenvalues()(i)>S) S=es.eigenvalues()(i);
				}
			RSO.Push(S);

			RSO.Push(S);
			o << 0 << "\t" << RSO.Mean()*3./2 << "\t" << RSO.StandardDeviation()*3./2 << "\n" ;
			RSO.Clear();


		for(int j=0;j<tMax;j+=Schrittweite){
			cout << "Schritt " << j << endl;

			if(MessePos==true) Positionen(ofbla,Chains);
			for(int k=0;k<Schrittweite*N*M;k++){

//				cout << "k=" << k << endl;
//                if(k%refreshrate==0) for(int j=0;j<M;j++) Chains->at(j).refreshNachbarn();

				Schritt=j*M*N+k;
				double phi=2*rndPI();
				double theta=acos(1-2*rnd());
				vektor<double> ev=vektor<double>(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta));
				ev=ev/ev.Betrag();
				int m=M*rnd();
				int i=N*rnd();
				Chains->at(m).EC(l,i,m,m,-1,ev);


				for(int i=0;i<M;i++)
					{
						vektor<double> q=Chains->at(i).Beads.at(0).ort-Chains->at(i).Beads.at(1).ort;
						q=q/q.Betrag();
						Q(0,0)+=q.x*q.x-1./3;
						Q(0,1)+=q.x*q.y;
						Q(0,2)+=q.x*q.z;
						Q(1,1)+=q.y*q.y-1./3;
						Q(1,2)+=q.y*q.z;
						Q(2,2)+=q.z*q.z-1./3;

				}
				Q(1,0)=Q(0,1);
				Q(2,0)=Q(0,2);
				Q(2,1)=Q(1,2);
				Q/=M;
		//		cout << "Q=" << Q << endl;
				SelfAdjointEigenSolver<Matrix3f> es(Q,false);
				double S=0;
				for(int i=0;i<3;i++){
						if(es.eigenvalues()(i)>S) S=es.eigenvalues()(i);
					}
				RSO.Push(S);
				}
			o << j+Schrittweite << "\t" << RSO.Mean()*3./2 << "\t" << RSO.StandardDeviation()*3./2 << "\n" ;
			RSO.Clear();
		}

	}

double Gesamtenergie(vector<chain> * Chains)
{
    double E=0;
    for(int i=0;i<Chains->size();i++)
    {
        for(int j=1;j<Chains->at(i).Beads.size();j++)
        {
            E+=Federenergie(Chains->at(i).Beads.at(j).ort-Chains->at(i).Beads.at(j-1).ort);
        }
    }

    return E/(M*(N-1));


}

void paarVerteilungBonds(vector<chain> *Chains, long long *rVal)
{
//    cout << "Funktion" << endl;
    int Genauigkeit=20;
    double maxAbstand=.5*sqrt(XMAX*XMAX+YMAX*YMAX+ZMAX*ZMAX);
    int Anz=0;
//    long long rVal[int(maxAbstand*Genauigkeit)];
//    for(int i=0;i<int(maxAbstand*Genauigkeit);i++) rVal[i]=0;
    for(int m=0;m<M;m++)
        {
        for(int n=0;n<N-1;n++)
        {
            vektor<double> b1=Chains->at(m).Beads.at(n+1).ort-Chains->at(m).Beads.at(n).ort;
            int j=n+1;
            for(int i=m+1;i<M;i++){

                for(j;j<N-1;j++)
                    {
//                        if(i==m) i++;
//                        if(i==M) break;
//                        if(i==m && j==n)
//                        {
//                            if(j==N-2){
//                                if(m==M-1) break;
//                                j=0;
//                                i++;
//                            }
//                            else j++;
//                        }
                        vektor<double> b2=Chains->at(i).Beads.at(j+1).ort-Chains->at(i).Beads.at(j).ort;
//                        int k0;
//                        int q0;
                        vektor<double> p1=Chains->at(m).Beads.at(n).ort;

                        double Abstand=1000;

                        for(int k=0;k<=20;k++){
//                            p2=Chains->at(i).Beads.at(j).ort;

                            double q=((Chains->at(m).Beads.at(n).ort-Chains->at(i).Beads.at(j).ort)*b2+b1*b2*1.*k/20)/(b2*b2);
                            if(q<0) q=0;
                            if(q>1) q=1;
                            vektor<double> p2=Chains->at(i).Beads.at(j).ort+b2*q;

//                            for(int q=0;q<=20;q++){


                                double s=(p1-p2).Betrag();
                                if(s<Abstand)  {
                                        Abstand=s;
//                                        k0=k;
//                                        q0=q;
                                }
                        //                                p2=p2+(b2*1./20);



                            p1=p1+(b1*1./20);

                        }
//                        vektor<double> b2=Chains->at(i).Beads.at(j+1).ort-Chains->at(i).Beads.at(j).ort;
////                        int k0;
////                        int q0;
//                        vektor<double> p1=Chains->at(m).Beads.at(n).ort;
//                        vektor<double> p2=Chains->at(i).Beads.at(j).ort;
//
//                        double Abstand=(p1-p2).Betrag();
//
//                        for(int k=0;k<=20;k++){
//                            p2=Chains->at(i).Beads.at(j).ort;
//
//                            for(int q=0;q<=20;q++){
//
//                                double s=(p1-p2).Betrag();
//                                if(s<Abstand)  {
//                                        Abstand=s;
////                                        k0=k;
////                                        q0=q;
//                                }
//                                p2=p2+(b2*1./20);
//
//
//                            }
//                            p1=p1+(b1*1./20);
//
//                        }
//                        cout << k0 << "\t" << q0 << endl;
//                        cout << "Abstand=" << Abstand << "\t erhöhe rVal an Stelle " <<  int(Abstand*Genauigkeit) << endl;
//                        if(int(Abstand*Genauigkeit)==0) cout << "Abstand=" << Abstand << "\t erhöhe rVal an Stelle " <<  int(Abstand*Genauigkeit) << endl;

//                        if(int(Abstand*Genauigkeit)==0) Anz++;
                        rVal[int(Abstand*Genauigkeit)]++;
                    }
                    j=0;
                }
        }
    }

//    cout << Anz << " Bonds in 0tem bin" << endl;
//    ofstream o;
//    o.open("paarVerteilungBonds.txt");
//    for(int i=0;i<int(XMAX/2*Genauigkeit);i++)
//        o << 1.*i/Genauigkeit+1./(2*Genauigkeit) << "\t" << 1.*rVal[i]/(N*N*M*M)*XMAX*YMAX*ZMAX/(4./3*PI*(pow(1.*(i+1)/Genauigkeit,3)-pow(1.*i/Genauigkeit,3))) << "\n" ;
//    o.close();

}
void paarVerteilungBeads(vector<chain> *Chains, long long *rVal)
{
//    cout << "Funktion" << endl;
    int Genauigkeit=20;
    double maxAbstand=.5*sqrt(XMAX*XMAX+YMAX*YMAX+ZMAX*ZMAX);
//    long long rVal[int(maxAbstand*Genauigkeit)];
    for(int m=0;m<M;m++)
        {
        for(int n=0;n<N;n++)
        {
            int j=n+1;
            for(int i=m+1;i<M;i++){

                for(j;j<N;j++)
                    {

                        double Abstand=(Chains->at(m).Beads.at(n).ort-Chains->at(i).Beads.at(j).ort).Betrag();
//                        cout << "Abstand=" << Abstand << endl;
                        rVal[int(Abstand*Genauigkeit)]++;
                    }
                    j=0;

                }
        }
    }

//    ofstream o;
//    o.open("paarVerteilungBeads.txt");
//    for(int i=0;i<int(XMAX/2*Genauigkeit);i++)
//        o << 1.*i/Genauigkeit+1./(2*Genauigkeit) << "\t" << 1.*rVal[i]/(N*N*M*M)*XMAX*YMAX*ZMAX/(4./3*PI*(pow(1.*(i+1)/Genauigkeit,3)-pow(1.*i/Genauigkeit,3))) << "\n" ;
//    o.close();

}

void heatmap(vector <chain> *Chains, int j1, int j2, double dist, ofstream &o)
{
    for(int i1=0;i1<N;i1++){
        for(int i2=0;i2<N;i2++){
            if((Chains->at(j1).Beads.at(i1).ort-Chains->at(j2).Beads.at(i2).ort).Betragsquadrat()<dist) o << Schritt/Schrittweite << "\t" << i1 << "\t" << i2 << endl;

        }
    }


}

void writeCheckpoint(vector<chain> *Chains)
{
    ofstream ofCheck;
    ofCheck.open("Checkpoint.txt");

    ofCheck.precision(17);

    ofCheck << Schritt << endl;
    ofCheck << AnzahlKollisionen << endl;
    ofCheck << AnzahlFederEC << endl;

    /**Positionen */
    for(int j=0;j<M;j++){
        for(int i =0;i<N;i++)
		    {
		       ofCheck << Chains->at(j).Beads.at(i).ort.x << "\t" << Chains->at(j).Beads.at(i).ort.y <<
	"\t"<< Chains->at(j).Beads.at(i).ort.z << endl;
	ofCheck << Chains->at(j).Beads.at(i).xtimes << "\t" << Chains->at(j).Beads.at(i).ytimes <<
	"\t"<< Chains->at(j).Beads.at(i).ztimes << endl;
	//               cout << k << " " << i << " " << j << endl;
		    }
//		    ofCheck << endl;
		}



    ofCheck << engine;

//    cout << Chains->at(0).midBead.size() << endl;

    cout << "Checkpoint fertig" << endl;
//    ofCheck.close();
}

vector<chain> * readCheckpoint(mt19937 *engine_,list<pair<int,int>> ***Zellen, string cp)
{
    vector<chain> *Chains = new vector<chain>();
    ifstream datei(cp,ios::in);
    if(!datei.is_open()) {
        cerr << "Datei " << cp << " konnte nicht geöffnet werden" << endl;
        throw exception();
    }
    double xVal[N*M];
    double yVal[N*M];
    double zVal[N*M];
    int xTimes[N*M];
    int yTimes[N*M];
    int zTimes[N*M];

    datei >> Schritt;
    datei >> AnzahlKollisionen;
    datei >> AnzahlFederEC;

    for(int i=0;i<N*M;i++)
    {
        datei >> xVal[i];
        datei >> yVal[i];
        datei >> zVal[i];
        datei >> xTimes[i];
        datei >> yTimes[i];
        datei >> zTimes[i];
    }



    for(int m=0;m<M;m++){

//                cout << midBead.size();

                Chains->push_back(chain(N,Chains,engine_,Zellen,&xVal[m*N],&yVal[m*N],&zVal[m*N],&xTimes[m*N],&yTimes[m*N],&zTimes[m*N]));

                Chains->at(m).engineX=engine_;
//                midBead.clear();

    }
    datei >> *engine_;





    return Chains;



}

