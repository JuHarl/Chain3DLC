#ifndef CHAIN_H
#define CHAIN_H


using namespace std;
struct bond{
    bead * Bead0;
    bead * Bead1;
    list<pair<int,int>>::iterator self;
    bond(bead * bead0, bead* bead1){
        Bead0=bead0;
        Bead1=bead1;
        Schwerpunkt=Bead0->ort+((Bead1->ort-Bead0->ort)*.5);
    }
    bond(){};
    vektor<double> Schwerpunkt;
    void refreshSchwerpunkt()
    {

        Schwerpunkt=Bead0->ort+((Bead1->ort-Bead0->ort)*.5);
        if(Schwerpunkt.x<0) Schwerpunkt.x+=XMAX;
        else if(Schwerpunkt.x>XMAX) Schwerpunkt.x-=XMAX;
        if(Schwerpunkt.y<0) Schwerpunkt.y+=YMAX;
        else if(Schwerpunkt.y>YMAX) Schwerpunkt.y-=YMAX;
        if(Schwerpunkt.z<0) Schwerpunkt.z+=ZMAX;
        else if(Schwerpunkt.z>ZMAX) Schwerpunkt.z-=ZMAX;
//        if(Schritt==Bug){
//            Bead0->ort.Ausgabe();
//            Bead1->ort.Ausgabe();
//            Schwerpunkt.Ausgabe();
//
//        }
    }
};

class chain
{
	public:
	int N;
	mt19937 *engineX;
    list<pair<int,int>> ***Zellen;
    vector<bond> Bonds;
	vector<bead> Beads;
	//vector<vektor<double>> Schwerpunkte;
	//vector<vektor<double>> midBead;
	//vector<vektor<double>> firstBead;
	//vector<vektor<double>> lastBead;
	vector<double> midBeadKorrMeans;
	vector<chain>* Chains;
	void NachbarDichte(double bin, int ChainIndex,ofstream &o);
	chain(int N, vector<chain>* Chains,mt19937 *engineX,list<pair<int,int>> ***Zellen);
	chain(int N, vector<chain>* Chains,mt19937 *engineX,list<pair<int,int>> ***Zellen, double * xVal, double * yVal, double * zVal,int *xTimes,int *yTimes,int *zTimes);
	double Gesamtenergie();
	//void Positionen(ofstream &ofbla);
	double ECMove(int i, int j, double x, vektor <double> ev);//Berechnet maximale Strecke die Bead i in Potential von Bead j zurücklegen kann
	void EC(double x,int i,int m,int lastchain, int lastbead, vektor<double> ev,ofstream &ofbla, bool beadkoll=false, bool rightKollBead=false);
	void EC(double x,int i,int m,int lastchain, int lastbead, vektor<double> ev,int d=0, bool beadkoll=false, bool rightKollBead=false, bool beadBeadKoll=false);
	void PosMidBead(ofstream &o);
	void KollisionsListe1(vector<pair<int,int>> * Koll1Bonds, int ** koll1,double x, int i,int m, vektor<double> ev);
	void KollisionsListe2(vector<pair<int,int>> * Koll2Bonds, int ** koll2,double x, int i,int m, vektor<double> ev);
	double maxWeg(vektor<double> i, vektor<double> j, vektor<double> ev,double deltaE);//Max Wegstrecke die bead i sich in Richtung ev bewegen darf, in Abhängigkeit der Energieänderung bzgl j
	double maxWegInnen(vektor<double> i, vektor<double> j, vektor<double> ev,double deltaE);
	double EndToEnd();
	//void Schwerpunkt(ofstream &o);
};


#endif //CHAIN_H
