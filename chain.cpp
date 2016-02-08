#include "chain.h"

chain::chain(int N, vector<chain>* Chains,mt19937 *engine ,  list<pair<int,int>> ***Zellen)
{
	engineX=engine;

    this->Zellen=Zellen;

	uniform_real_distribution<> distribution(0, 1);
	auto rnd = bind(distribution, ref(*engineX));

	uniform_real_distribution<> PIdistribution(0, PI);
	auto rndPI = bind(PIdistribution, ref(*engineX));

	uniform_real_distribution<> NEGdistribution(-1, 1);
	auto rndWAY = bind(NEGdistribution, ref(*engineX));


//	cout << "Konstruktor Chain" << endl;
	this->Chains=Chains;
//	cout << "Pushe in Chains" << endl;
//	Chains->push_back(*this);
	this->N=N;
	vektor<double> x;
	x=vektor<double>(rnd()*XMAX,rnd()*YMAX,rnd()*ZMAX);
//	x=vektor<double>(.5*XMAX,.5*YMAX,.5*ZMAX);
//	if(k==0) x=vektor<double>(.5*XMAX+1.5*x0,.5*YMAX+x0/2+.1,.5*ZMAX-1.*x0);
//	else x=vektor<double>(.5*XMAX+x0,.5*YMAX-x0/2,.5*ZMAX-.5*x0);
//	cout << "..." << endl;
	Beads.push_back(x);

    vektor<double> ex;
//    int xtimes=0;
//    int ytimes=0;
//    int ztimes=0;
	for(int i=1;i<N;i++)
		{
//			cout << "Bead " << i << endl;
			double phi=2*rndPI();
			double theta=acos(1-2*rnd());
			ex=vektor<double>(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
//			ex=vektor<double>(sin(phi),cos(phi),1);
//			if(k==0) ex=vektor<double>(-1,0,0);
//			else ex=vektor<double>(0,0,-1);
//			ex=vektor<double>(1,0,0);
//			if(i<N/4) ex=vektor<double>(-1,0,0);
//			if(i==N/4) ex=vektor<double>(0,-1,0);
//			if(i>N/4 && i< N/2) ex=vektor<double>(1,0,1);
//			if(i>=N/2 && i<3*N/4) ex=vektor<double>(0,0,-1);
//			if(i==3*N/4 || i==3*N/4+1) ex=vektor<double>(0,1,0);
//			if(i>3*N/4+1) ex=vektor<double>(0,0,1);

			ex=ex/ex.Betrag()*x0;



			//vektor<double> x=Beads.at(i-1).ort-ex;
			Beads.push_back(bead(Beads.at(i-1).ort));
//			vektor<double> newOrt=Beads.at(i-1).ort.add(ex);
//			if (newOrt.x>XMAX) xtimes++;
//			else if(newOrt.x<0) xtimes--;
//			if (newOrt.y>YMAX) ytimes++;
//			else if(newOrt.y<0) ytimes--;
//			if (newOrt.x>XMAX) xtimes++;
//			else if(newOrt.x<0) xtimes--;
			Beads.at(i).verschiebe(ex);
			Beads.at(i).xtimes+=Beads.at(i-1).xtimes;
			Beads.at(i).ytimes+=Beads.at(i-1).ytimes;
			Beads.at(i).ztimes+=Beads.at(i-1).ztimes;
		}

    for(int i=0;i<N-1;i++)
    {
        Bonds.push_back(bond(&Beads.at(i),&Beads.at(i+1)));
        Zellen[int(Bonds.at(i).Schwerpunkt.x/(Zellgroesse))][int(Bonds.at(i).Schwerpunkt.y/(Zellgroesse))][int(Bonds.at(i).Schwerpunkt.z/(Zellgroesse))].push_back(make_pair(Chains->size(),i));
        Bonds.at(i).self=Zellen[int(Bonds.at(i).Schwerpunkt.x/(Zellgroesse))][int(Bonds.at(i).Schwerpunkt.y/(Zellgroesse))][int(Bonds.at(i).Schwerpunkt.z/(Zellgroesse))].end();
        Bonds.at(i).self--;
//        cout << "Füge Schwerpunkt";
//        Bonds.at(i).Schwerpunkt.Ausgabe();
//        cout << "in Zelle" << int(Bonds.at(i).Schwerpunkt.x)/(Zellgroesse) << int(Bonds.at(i).Schwerpunkt.y)/(Zellgroesse) <<int(Bonds.at(i).Schwerpunkt.z)/(Zellgroesse) << "ein" << endl;
    }
//    cout << "maxSize=" <<  midBead.max_size() << endl;




}

chain::chain(int N, vector<chain>* Chains,mt19937 *engine,list<pair<int,int>> ***Zellen, double * xVal, double * yVal, double * zVal,int *xTimes,int *yTimes,int *zTimes)
{
	engineX=engine;

    this->Zellen=Zellen;

	uniform_real_distribution<> distribution(0, 1);
	auto rnd = bind(distribution, ref(*engineX));

	uniform_real_distribution<> PIdistribution(0, PI);
	auto rndPI = bind(PIdistribution, ref(*engineX));

	uniform_real_distribution<> NEGdistribution(-1, 1);
	auto rndWAY = bind(NEGdistribution, ref(*engineX));

	uniform_int_distribution<> INTdistribution(0, N-1);
	auto INTrnd = bind(INTdistribution, ref(*engineX));
//	cout << "Konstruktor Chain" << endl;
	this->Chains=Chains;
//	cout << "Pushe in Chains" << endl;
//	Chains->push_back(*this);
	this->N=N;
	vektor<double> x;
//	x=vektor<double>(rnd()*XMAX,rnd()*YMAX,rnd()*ZMAX);
//	x=vektor<double>(.5*XMAX,.5*YMAX,.5*ZMAX);
//	if(k==0) x=vektor<double>(.5*XMAX+1.5*x0,.5*YMAX+x0/2+.1,.5*ZMAX-1.*x0);
//	else x=vektor<double>(.5*XMAX+x0,.5*YMAX-x0/2,.5*ZMAX-.5*x0);
//	cout << "..." << endl;
//	Beads.push_back(x);

//    vektor<double> ex;

	for(int i=0;i<N;i++)
		{
//			cout << "Bead " << i << endl;
//			double phi=2*rndPI();
//			double theta=acos(1-2*rnd());
//			ex=vektor<double>(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
//			ex=vektor<double>(sin(phi),cos(phi),1);
//			if(k==0) ex=vektor<double>(-1,0,0);
//			else ex=vektor<double>(0,0,-1);
//			ex=vektor<double>(1,0,0);
//			if(i<N/4) ex=vektor<double>(-1,0,0);
//			if(i==N/4) ex=vektor<double>(0,-1,0);
//			if(i>N/4 && i< N/2) ex=vektor<double>(1,0,1);
//			if(i>=N/2 && i<3*N/4) ex=vektor<double>(0,0,-1);
//			if(i==3*N/4 || i==3*N/4+1) ex=vektor<double>(0,1,0);
//			if(i>3*N/4+1) ex=vektor<double>(0,0,1);
//			ex=ex/ex.Betrag()*x0;
			//vektor<double> x=Beads.at(i-1).ort-ex;
			Beads.push_back(bead(xVal[i],yVal[i],zVal[i],xTimes[i],yTimes[i],zTimes[i]));
//			Beads.at(i).verschiebe(ex);
		}

    for(int i=0;i<N-1;i++)
    {
        Bonds.push_back(bond(&Beads.at(i),&Beads.at(i+1)));
        Zellen[int(Bonds.at(i).Schwerpunkt.x/(Zellgroesse))][int(Bonds.at(i).Schwerpunkt.y/(Zellgroesse))][int(Bonds.at(i).Schwerpunkt.z/(Zellgroesse))].push_back(make_pair(Chains->size(),i));
        Bonds.at(i).self=Zellen[int(Bonds.at(i).Schwerpunkt.x/(Zellgroesse))][int(Bonds.at(i).Schwerpunkt.y/(Zellgroesse))][int(Bonds.at(i).Schwerpunkt.z/(Zellgroesse))].end();
        Bonds.at(i).self--;
//        cout << "Füge Schwerpunkt";
//        Bonds.at(i).Schwerpunkt.Ausgabe();
//        cout << "in Zelle" << int(Bonds.at(i).Schwerpunkt.x)/(Zellgroesse) << int(Bonds.at(i).Schwerpunkt.y)/(Zellgroesse) <<int(Bonds.at(i).Schwerpunkt.z)/(Zellgroesse) << "ein" << endl;
    }




}

void chain::NachbarDichte(double bin, int ChainIndex,ofstream &o)
{
        int Anz[int(XMAX/bin)][int(YMAX/bin)][int(ZMAX/bin)];
        for(int x=0;x<int(XMAX/bin);x++)
        for(int y=0;y<int(YMAX/bin);y++)
        for(int z=0;z<int(ZMAX/bin);z++)
            Anz[x][y][z]=0;

        for(int i=0;i<N;i++){
        int Xmin=(int(Beads.at(i).ort.x/Zellgroesse)-1+AnzZellen)%AnzZellen;
        int Xmax=(int(Beads.at(i).ort.x/Zellgroesse)+2)%AnzZellen;
        int Ymin=(int(Beads.at(i).ort.y/Zellgroesse)-1+AnzZellen)%AnzZellen;
        int Ymax=(int(Beads.at(i).ort.y/Zellgroesse)+2)%AnzZellen;
        int Zmin=(int(Beads.at(i).ort.z/Zellgroesse)-1+AnzZellen)%AnzZellen;
        int Zmax=(int(Beads.at(i).ort.z/Zellgroesse)+2)%AnzZellen;


                for(int X=Xmin;X!=Xmax;X=(X+1)%AnzZellen){
                for(int Y=Ymin;Y!=Ymax;Y=(Y+1)%AnzZellen){
                for(int Z=Zmin;Z!=Zmax;Z=(Z+1)%AnzZellen){
                        for(list<pair<int,int>>::iterator it=Zellen[X][Y][Z].begin();it!=Zellen[X][Y][Z].end();++it)
                            {
                                int j=it->first;
                                int k=it->second;
                                if(j!=ChainIndex){
                                    if((Chains->at(j).Beads.at(k).ort-Beads.at(i).ort).Betragsquadrat()<bin)
                                    Anz[int(Chains->at(j).Beads.at(k).ort.x/bin)][int(Chains->at(j).Beads.at(k).ort.y/bin)][int(Chains->at(j).Beads.at(k).ort.z/bin)]++;
//                                    o << Chains->at(j).Beads.at(k).ort.x << "\t" << Chains->at(j).Beads.at(k).ort.y << "\t" << Chains->at(j).Beads.at(k).ort.z << "\n";
                                }
                            }
                }
                }
                }
        }
    for(int x=0;x<int(XMAX/bin);x++)
    for(int y=0;y<int(YMAX/bin);y++)
    for(int z=0;z<int(ZMAX/bin);z++)
        if(Anz[x][y][z]!=0) o << 1.*x/bin << "\t" << 1.*y/bin << "\t" << 1.*z/bin << "\t" << Anz[x][y][z] << endl;



}
double chain::EndToEnd()
	{
	    vektor<double> addRB=vektor<double>((Beads.at(N-1).xtimes-Beads.at(0).xtimes)*XMAX,(Beads.at(N-1).ytimes-Beads.at(0).ytimes)*YMAX,(Beads.at(N-1).ztimes-Beads.at(0).ztimes)*ZMAX);
		return (Beads.at(N-1).ort.sub(Beads.at(0).ort).add(addRB)).Betragsquadrat();
	}

void chain::PosMidBead(ofstream &o)//Speichert die Position des mittleren beads
{
	//cout << "PushmidBead" << endl;
	//vektor<double> x=vektor<double>(Beads.at(N/2).ort.x,Beads.at(N/2).ort.y,Beads.at(N/2).ort.z);
    o.precision(100);

	vektor<double> S= vektor<double>(0,0,0);
		for(int i=0;i<N;i++)
			{
				S.x+=Beads.at(i).ort.x+Beads.at(i).xtimes*XMAX;
				S.y+=Beads.at(i).ort.y+Beads.at(i).ytimes*YMAX;
				S.z+=Beads.at(i).ort.z+Beads.at(i).ztimes*ZMAX;
			}
    o << (S.x/N) << "\t"<< (S.y/N) << "\t"<< (S.z/N) << "\n";

	vektor<double> x=vektor<double>(Beads.at(N/2).ort.x+Beads.at(N/2).xtimes*XMAX,Beads.at(N/2).ort.y+Beads.at(N/2).ytimes*YMAX,Beads.at(N/2).ort.z+Beads.at(N/2).ztimes*ZMAX);
	//midBead.push_back(x);
	o << x.x << "\t"<< x.y << "\t"<< x.z << "\n";
    x=vektor<double>(Beads.at(0).ort.x+Beads.at(0).xtimes*XMAX,Beads.at(0).ort.y+Beads.at(0).ytimes*YMAX,Beads.at(0).ort.z+Beads.at(0).ztimes*ZMAX);
	//firstBead.push_back(x);
	o << x.x << "\t"<< x.y << "\t"<< x.z<< "\n";
	x=vektor<double>(Beads.at(N-1).ort.x+Beads.at(N-1).xtimes*XMAX,Beads.at(N-1).ort.y+Beads.at(N-1).ytimes*YMAX,Beads.at(N-1).ort.z+Beads.at(N-1).ztimes*ZMAX);
	//lastBead.push_back(x);
	o << x.x << "\t"<< x.y << "\t"<< x.z<< "\n";



}


double chain::Gesamtenergie()
{
	double E=0;

	for(int i=0;i<N-1;i++){
		E+=K/2*pow((Beads.at(i).ort-Beads.at(i+1).ort).Betrag()-x0,2);
	}
	return E/N;
}



void chain::EC(double x, int i,int m,int lastchain, int lastbead, vektor<double> ev, int d, bool beadKoll, bool rightKollBead, bool beadBeadKoll)
{
//    cout.precision(20);

//	uniform_real_distribution<> distribution(0, 1);
//	auto rnd = bind(distribution, ref(*engineX));
//
//	uniform_real_distribution<> PIdistribution(0, PI);
//	auto rndPI = bind(PIdistribution, ref(*engineX));
//
//	uniform_real_distribution<> NEGdistribution(-1, 1);
//	auto rndWAY = bind(NEGdistribution, ref(*engineX));
//
//	uniform_int_distribution<> INTdistribution(0, N-1);
//	auto INTrnd = bind(INTdistribution, ref(*engineX));
//
//	uniform_int_distribution<> Mdistribution(0, M-1);
//	auto Mrnd = bind(Mdistribution, ref(*engineX));
//    if(d==0) cout << "RANDOMZAHL IN EC=" << rnd() << endl;

	if(x<=0){
		if(x<0) cout << "----Fehler: x=" << x << "!!!--- \n EC beendet" << endl;
//		if(Schritt==Bug) cout << "EC beendet" << endl;
		return;
	}
	if(d>1e4){//Einfache Lösung, noch zu überprüfen
		cout << "Endlosschleife gefunden, fahre mit EC fort, x=" << x << endl;
		AnzahlEndlosschleifen++;
        double phi=2*rndPI();
        double theta=acos(1-2*rnd());
        vektor<double> temp=vektor<double>(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta));
        ev=ev.kreuz(temp);//Randomvektor, senkrecht zu ev.
        ev=ev/ev.Betrag();
        d=0;
//		int nextChain=M*rnd();
//		Chains->at(nextChain).EC(x,N*rnd(),nextChain,-1,-1,ev,0);
//		return;
	}
//	if(Schritt==Bug) cout << "\n \n \n Bewegung von Kette " << m << " Bead " << i << " x=" << fixed << x << " d=" << d << endl;
//	if(Schritt==Bug) cout << "lastchain =" << lastchain << " lastbead=" << lastbead << endl;
//	if(Schritt==Bug && beadKoll) cout << "lastbead= " << lastbead+1 << endl;
//	if(Schritt==Bug && rightKollBead) cout << "rightKollBead " << endl;
//	if(Schritt==Bug && beadBeadKoll) cout << "beadBeadKoll " << endl;
//	double phi=2*rndPI();
//	double theta=rndPI();
	//vektor<double> ev=vektor<double>(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta));
	//vektor<double> ev=vektor<double>(0,0,1);
	//ev=ev/ev.Betrag();
	//ev.Ausgabe();
	//if(ev.Betrag()!=1) {cout << "Fehler: Betrag von ev ist " << ev.Betrag() << endl;}
	//double x=l;

	//while(x>0){
	double w=x;
	double q=x;

	bool thisBeadKoll;
	bool nextBeadKoll=false;

	int nextKollChain;
	int nextKollBead;
	bool rightchain;

//
//    cout << "Erstelle Kollisionsliste" << endl;



//    vector<pair<int,int>> Koll1Bonds;
//    vector<pair<int,int>> Koll2Bonds;

//    if(Schritt==Bug)
//    {
//        ListeAusgeben(listX,0);
//        ListeAusgeben(listY,1);
//        ListeAusgeben(listZ,2);
//
//        cout << "Bewegung von ";
//        Beads.at(i).ort.Ausgabe();
//        cout << " nach";
//        (Beads.at(i).ort.add(ev*x)).Ausgabe();
//        cout << "ev=";
//        ev.Ausgabe();
//        cout << endl;
//    }
//    int mglKoll=0;
//    int tatKoll=0;



//    cout << "Berechne Kollisionen" << endl;



//					if(Schritt==Bug) cout << "\n Prüfe auf Kollision mit Kette " << j <<  " Beads " << k  << "," << k-1 << endl;
					double q1=x;
					double q2=x;
					double p1=1000;
					double p2=1000;
					double p=1000;
					int j1,k1,j2,k2;
					bool nextBeadKoll1=false;
					bool nextBeadKoll2=false;
//					if (Schritt==Bug) cout << "linker Bond" << endl;
					int xSchwerL;
					int ySchwerL;
					int zSchwerL;

					if(i!=0){
//                        if(Schritt==Bug){
//                        cout << "Schwerpunkt" ;
//                        Bonds.at(i-1).Schwerpunkt.Ausgabe();
//                        cout << "X=" << ((int(Bonds.at(i-1).Schwerpunkt.x)/Zellgroesse)-1+AnzZellen)%AnzZellen << endl;
//                        cout << "XMax=" << (int(Bonds.at(i-1).Schwerpunkt.x)/Zellgroesse+2)%AnzZellen << endl;
//                        cout << "Y=" << ((int(Bonds.at(i-1).Schwerpunkt.y)/Zellgroesse)-1+AnzZellen)%AnzZellen << endl;
//                        cout << "YMax=" << (int(Bonds.at(i-1).Schwerpunkt.y)/Zellgroesse+2)%AnzZellen << endl;
//                        cout << "Z=" << ((int(Bonds.at(i-1).Schwerpunkt.z)/Zellgroesse)-1+AnzZellen)%AnzZellen << endl;
//                        cout << "ZMax=" << (int(Bonds.at(i-1).Schwerpunkt.z)/Zellgroesse+2)%AnzZellen << endl;
                        vektor<double> r_=Beads.at(i).ort-Beads.at(i-1).ort;
                        xSchwerL=Bonds.at(i-1).Schwerpunkt.x/Zellgroesse;
                        ySchwerL=Bonds.at(i-1).Schwerpunkt.y/Zellgroesse;
                        zSchwerL=Bonds.at(i-1).Schwerpunkt.z/Zellgroesse;
                        int Xmin=(xSchwerL-1+AnzZellen)%AnzZellen;
                        int Xmax=(xSchwerL+2)%AnzZellen;
                        int Ymin=(ySchwerL-1+AnzZellen)%AnzZellen;
                        int Ymax=(ySchwerL+2)%AnzZellen;
                        int Zmin=(zSchwerL-1+AnzZellen)%AnzZellen;
                        int Zmax=(zSchwerL+2)%AnzZellen;

//                        if(x>x0){
//                            if(ev.x>0) Xmax=(Xmax+1)%AnzZellen;
//                            else Xmin=(Xmin-1+AnzZellen)%AnzZellen;
//                            if(ev.y>0) Ymax=(Ymax+1)%AnzZellen;
//                            else Ymin=(Ymin-1+AnzZellen)%AnzZellen;
//                            if(ev.z>0) Zmax=(Zmax+1)%AnzZellen;
//                            else Zmin=(Zmin-1+AnzZellen)%AnzZellen;
//                        }
//
                        for(int X=Xmin;X!=Xmax;X=(X+1)%AnzZellen){
                        for(int Y=Ymin;Y!=Ymax;Y=(Y+1)%AnzZellen){
                        for(int Z=Zmin;Z!=Zmax;Z=(Z+1)%AnzZellen){
//                                cout <<int(Bonds.at(i-1).Schwerpunkt.y)%Zellgroesse << endl;
//                                if(Schritt==Bug) cout << X << Y << Z << endl;
//                        cout << Zellen[X][Y][Z].size() << endl;
//                        if(!Zellen[X][Y][Z].empty()) {
                            for(list<pair<int,int>>::iterator it=Zellen[X][Y][Z].begin();it!=Zellen[X][Y][Z].end();++it)
                            {

//                                if(Schritt==Bug) cout << "for1 ";

                                int j=it->first;
                                int k=it->second;
//                                if(Schritt==Bug) cout << "j="<<j << "k=" << k << endl;
                                if((j==m && (k==i-2 || k==i || k==i-1)) || (j==lastchain && k==lastbead-1 && (rightKollBead || beadBeadKoll)) || (j==lastchain && k==lastbead && (rightKollBead || beadBeadKoll) && beadKoll) );
                                else {
            //					    if(Schritt==Bug) cout << "Prüfe ausgehend von Beads "<< i-1 << "," << i << endl;
                                    pair<double,double> Kollp=Kollision(Beads.at(i-1).ort, Chains->at(j).Bonds.at(k).Bead0->ort, Chains->at(j).Bonds.at(k).Bead1->ort, r_, ev,x); //Überprüfung des linken Bonds
//                                    if(Schritt==Bug) cout << "w=" << Koll << endl;
//                                    mglKoll++;

//                                    if(Koll<x) tatKoll++;
                                    double Koll=Kollp.first;


                                    if(abs(Koll-q1)<1e-10 && abs(Koll-x)>1e-10 && abs(k-k1)<=1 && j==j1)//Bead wird direkt getroffen
                                    {
        //        						if(Schritt==Bug) cout << "Bead wird direkt getroffen" << endl;
                                        nextBeadKoll1=true;
                                        k=min(k,k1);
                                        k1=k;
        //                                if(Schritt==Bug) cout << "k=" << k << endl;
                                    }
                                    if(Koll < q1) {
                                            q1=Koll;
                                            p1=Kollp.second;
                                            j1=j;
                                            k1=k;
                                            if(abs(Koll-q1)>1e-10) nextBeadKoll1=false;

                                            }

                            }
                        }
                        }}}
//                        }
					}
//					cout << "Rechter Bond " << endl;
                    int xSchwerR;
                    int ySchwerR;
                    int zSchwerR;
					if(i!=N-1){
                        vektor<double> r_=Beads.at(i).ort-Beads.at(i+1).ort;


                        xSchwerR=Bonds.at(i).Schwerpunkt.x/Zellgroesse;
                        ySchwerR=Bonds.at(i).Schwerpunkt.y/Zellgroesse;
                        zSchwerR=Bonds.at(i).Schwerpunkt.z/Zellgroesse;
                        int Xmin=(xSchwerR-1+AnzZellen)%AnzZellen;
                        int Xmax=(xSchwerR+2)%AnzZellen;
                        int Ymin=(ySchwerR-1+AnzZellen)%AnzZellen;
                        int Ymax=(ySchwerR+2)%AnzZellen;
                        int Zmin=(zSchwerR-1+AnzZellen)%AnzZellen;
                        int Zmax=(zSchwerR+2)%AnzZellen;

//                        if(x>x0){
//                            if(ev.x>0) Xmax=(Xmax+1)%AnzZellen;
//                            else Xmin=(Xmin-1+AnzZellen)%AnzZellen;
//                            if(ev.y>0) Ymax=(Ymax+1)%AnzZellen;
//                            else Ymin=(Ymin-1+AnzZellen)%AnzZellen;
//                            if(ev.z>0) Zmax=(Zmax+1)%AnzZellen;
//                            else Zmin=(Zmin-1+AnzZellen)%AnzZellen;
//                        }
////
//                        }
                        for(int X=Xmin;X!=Xmax;X=(X+1)%AnzZellen){
                        for(int Y=Ymin;Y!=Ymax;Y=(Y+1)%AnzZellen){
                        for(int Z=Zmin;Z!=Zmax;Z=(Z+1)%AnzZellen){
//                            if(Schritt==Bug) cout << X << Y << Z << endl;

//                        if(!Zellen[X][Y][Z].empty()){
                                for(list<pair<int,int>>::iterator it=Zellen[X][Y][Z].begin();it!=Zellen[X][Y][Z].end();++it)                        {
//                            if(Schritt==Bug) cout << "for2 ";

                            int j=it->first;
                            int k=it->second;
//                            if(Schritt==Bug) cout << "j="<<j << "k=" << k << endl;

                            if((j==m && (k==i+1 ||k==i || k==i-1)) || (j==lastchain && k==lastbead-1 && (!rightKollBead|| beadBeadKoll)) || (j==lastchain && k==lastbead && (!rightKollBead || beadBeadKoll) && beadKoll) );
                            else {

        //					    if(Schritt==Bug) cout << "Prüfe ausgehend von Beads "<< i-1 << "," << i << endl;
                                pair<double,double> Kollp=Kollision(Beads.at(i+1).ort, Chains->at(j).Bonds.at(k).Bead0->ort, Chains->at(j).Bonds.at(k).Bead1->ort, r_, ev,x); //Überprüfung des linken Bonds
//                                if(Schritt==Bug) cout << "w=" << Koll << endl;
//                                mglKoll++;
                                double Koll=Kollp.first;
//                                if(Koll<x) tatKoll++;
                                if(abs(Koll-q2)<1e-10 && abs(Koll-x)>1e-10 && abs(k-k2)<=1 && j==j2)//Bead wird direkt getroffen
                                {
    //        						if(Schritt==Bug) cout << "Bead wird direkt getroffen" << endl;
                                    nextBeadKoll2=true;
                                    k=min(k,k2);
                                    k2=k;
                                }
                                if(Koll < q2) {
                                        q2=Koll;
                                        p2=Kollp.second;
                                        j2=j;
                                        k2=k;
                                        if(abs(Koll-q2)>1e-10) nextBeadKoll2=false;

                                    }
                            }
                        } }}}
//                        }
					}
//					if(Schritt==Bug) cout << "q1=" << q1 << " j1=" << j1 << " k1=" << k1<< " q2=" << q2 << " j2=" << j2 << " k2=" << k2 << endl;
					q=min(q1,q2);
					//cout << "---" << endl;
					if(q<0) cout << "--------------FEHLER, q<0! q=" << q << "---------------" << endl;


					if(q<w) //Kollision mit Chain j Beads k-1 & k
						{
						    int j,k;
							if(abs(q1-q2)<1e-10 && x-q1>1e-10 &&j1==j2 && abs(k1-k2)<=1) {
//								if(Schritt==Bug) cout << "Kollision direkt mit Bead!" << endl;
								thisBeadKoll=true; //Bond wird direkt von bead getroffen
								j=j1;
								k=min(k1,k2);
								nextBeadKoll=(nextBeadKoll1 || nextBeadKoll2);
							}
							else {thisBeadKoll=false;
                            if(q2<q1) {
                                nextBeadKoll=nextBeadKoll2;
                                rightchain=true;
                                j=j2;
                                k=k2;
                                p=p2;
							}
							else {
							    nextBeadKoll=nextBeadKoll1;
                                rightchain=false;
                                j=j1;
                                k=k1;
                                p=p1;
							}
							}
//                            if(Schritt==Bug) cout << "k=" << k << endl;

							w=q;
							nextKollChain=j;
							nextKollBead=k+1;
//							nextBeadKoll=false;
						}




//	if(w!=x) w-=1e-3;// Bewegung nicht ganz bis zum Bond! ACHTUNG!
//	cout << "x= " << x << endl;
	double w1=ECMove(i,i+1,x,ev);
	double w2=ECMove(i,i-1,x,ev);



    double wMin=min(w,min(w1,w2));

//    cout << "Anteil Kollision=" << 1.*tatKoll/mglKoll << endl;





//    cout << "verschiebe" << endl;
    Beads.at(i).verschiebe(ev*wMin);
//    if(Schritt==Bug) cout << "aktualisiere linken Bond" << endl;
    if(i!=0) {
//            if(Schritt==Bug) cout << "Setze Zeiger auf alte Zelle" << endl;


//            int Xold=int(Bonds.at(i-1).Schwerpunkt.x/(Zellgroesse));
//            int Yold=int(Bonds.at(i-1).Schwerpunkt.y/(Zellgroesse));
//            int Zold=int(Bonds.at(i-1).Schwerpunkt.z/(Zellgroesse));

            Bonds.at(i-1).refreshSchwerpunkt();


//            if(Schritt==Bug) Bonds.at(i).Schwerpunkt.Ausgabe();
//            if(Schritt==Bug) cout << "neue Zelle " << int(Bonds.at(i).Schwerpunkt.x)/(Zellgroesse) << int(Bonds.at(i).Schwerpunkt.y)/(Zellgroesse) << int(Bonds.at(i).Schwerpunkt.z)/(Zellgroesse) << endl;
//            if(Schritt==Bug) cout << " hat " << Zellen[int(Bonds.at(i).Schwerpunkt.x)/(Zellgroesse)][int(Bonds.at(i).Schwerpunkt.y)/(Zellgroesse)][int(Bonds.at(i).Schwerpunkt.z)/(Zellgroesse)].size() << " Elemente" << endl;

            int X=int(Bonds.at(i-1).Schwerpunkt.x/(Zellgroesse));
            int Y=int(Bonds.at(i-1).Schwerpunkt.y/(Zellgroesse));
            int Z=int(Bonds.at(i-1).Schwerpunkt.z/(Zellgroesse));
            if(X>AnzZellen-1 || X<0) cout << "Fehler! Falsche Zellenzuordnung. X=" << X << endl;
            if(Y>AnzZellen-1 || Y<0) cout << "Fehler! Falsche Zellenzuordnung. X=" << Y << endl;
            if(Z>AnzZellen-1 || Z<0) cout << "Fehler! Falsche Zellenzuordnung. X=" << Z << endl;

            if(X!=xSchwerL || Y!=ySchwerL || Z!=zSchwerL){
                list<pair<int,int>> * old = &Zellen[xSchwerL][ySchwerL][zSchwerL];


//            if(Schritt==Bug) Bonds.at(i-1).Schwerpunkt.Ausgabe();
//            if(Schritt==Bug) cout << "neue Zelle " << int(Bonds.at(i-1).Schwerpunkt.x)/(Zellgroesse) << int(Bonds.at(i-1).Schwerpunkt.y)/(Zellgroesse) << int(Bonds.at(i-1).Schwerpunkt.z)/(Zellgroesse) << endl;
//            if(Schritt==Bug) cout << " hat " << Zellen[int(Bonds.at(i-1).Schwerpunkt.x)/(Zellgroesse)][int(Bonds.at(i-1).Schwerpunkt.y)/(Zellgroesse)][int(Bonds.at(i-1).Schwerpunkt.z)/(Zellgroesse)].size() << " Elemente" << endl;

                list<pair<int,int>>::iterator einf = Zellen[X][Y][Z].begin();
    //            cout  << Bonds.at(i-1).self->second << endl;
    //            if(Schritt==Bug) cout << "splice" << endl;

                Zellen[X][Y][Z].splice(einf,*old,Bonds.at(i-1).self);

    //            if(Schritt==Bug) cout << "Setze neuen iterator" << endl;

                Bonds.at(i-1).self=Zellen[X][Y][Z].begin();
//                Bonds.at(i-1).self;
            }

    }
//    if(Schritt==Bug) cout << "aktualisiere rechten Bond" << endl;
    if(i!=N-1){
//            if(Schritt==Bug) cout << "Setze Zeiger auf alte Zelle" << endl;


//            if(Schritt==Bug) cout << "alte Zelle " << int(Bonds.at(i).Schwerpunkt.x)/(Zellgroesse) << int(Bonds.at(i).Schwerpunkt.y)/(Zellgroesse) << int(Bonds.at(i).Schwerpunkt.z)/(Zellgroesse) << endl;
//            if(Schritt==Bug) Bonds.at(i).Schwerpunkt.Ausgabe();
//            int Xold=int(Bonds.at(i).Schwerpunkt.x/(Zellgroesse));
//            int Yold=int(Bonds.at(i).Schwerpunkt.y/(Zellgroesse));
//            int Zold=int(Bonds.at(i).Schwerpunkt.z/(Zellgroesse));

            Bonds.at(i).refreshSchwerpunkt();


//            if(Schritt==Bug) Bonds.at(i).Schwerpunkt.Ausgabe();
//            if(Schritt==Bug) cout << "neue Zelle " << int(Bonds.at(i).Schwerpunkt.x)/(Zellgroesse) << int(Bonds.at(i).Schwerpunkt.y)/(Zellgroesse) << int(Bonds.at(i).Schwerpunkt.z)/(Zellgroesse) << endl;
//            if(Schritt==Bug) cout << " hat " << Zellen[int(Bonds.at(i).Schwerpunkt.x)/(Zellgroesse)][int(Bonds.at(i).Schwerpunkt.y)/(Zellgroesse)][int(Bonds.at(i).Schwerpunkt.z)/(Zellgroesse)].size() << " Elemente" << endl;

            int X=int(Bonds.at(i).Schwerpunkt.x/(Zellgroesse));
            int Y=int(Bonds.at(i).Schwerpunkt.y/(Zellgroesse));
            int Z=int(Bonds.at(i).Schwerpunkt.z/(Zellgroesse));
            if(X>AnzZellen-1 || X<0) cout << "Fehler! Falsche Zellenzuordnung. X=" << X << endl;
            if(Y>AnzZellen-1 || Y<0) cout << "Fehler! Falsche Zellenzuordnung. X=" << Y << endl;
            if(Z>AnzZellen-1 || Z<0) cout << "Fehler! Falsche Zellenzuordnung. X=" << Z << endl;

            if(X!=xSchwerR || Y!=ySchwerR || Z!=zSchwerR){
                list<pair<int,int>> * old = &Zellen[xSchwerR][ySchwerR][zSchwerR];


//            if(Schritt==Bug) Bonds.at(i-1).Schwerpunkt.Ausgabe();
//            if(Schritt==Bug) cout << "neue Zelle " << int(Bonds.at(i-1).Schwerpunkt.x)/(Zellgroesse) << int(Bonds.at(i-1).Schwerpunkt.y)/(Zellgroesse) << int(Bonds.at(i-1).Schwerpunkt.z)/(Zellgroesse) << endl;
//            if(Schritt==Bug) cout << " hat " << Zellen[int(Bonds.at(i-1).Schwerpunkt.x)/(Zellgroesse)][int(Bonds.at(i-1).Schwerpunkt.y)/(Zellgroesse)][int(Bonds.at(i-1).Schwerpunkt.z)/(Zellgroesse)].size() << " Elemente" << endl;

                list<pair<int,int>>::iterator einf = Zellen[X][Y][Z].begin();
    //            cout  << Bonds.at(i-1).self->second << endl;
    //            if(Schritt==Bug) cout << "splice" << endl;

                Zellen[X][Y][Z].splice(einf,*old,Bonds.at(i).self);

    //            if(Schritt==Bug) cout << "Setze neuen iterator" << endl;

                Bonds.at(i).self=Zellen[X][Y][Z].begin();
//                Bonds.at(i).self;
            }
    }

//	if(Schritt==Bug) cout << "===w1= " << w1 << " ; " << " w2= " << w2 << " wKoll="<< fixed  << w << "===" << endl;
	if(w<w1 && w < w2){
//            cout << "Kollision" << endl;

//			if(Schritt==Bug) cout << "Kollision Chain " << m << " mit " << nextKollChain << endl;
			AnzahlKollisionen++;
//			if(Schritt==Bug) cout << "Bead " << i << " wird um w" << w << "bewegt" << endl;
//			Beads.at(i).verschiebe(ev*w);
//			if(Schritt==Bug){
//                if(rightchain) cout << "Kollision ausgehend von rechtem Bond" << endl;
//                else cout << "Kollision ausgehend von linkem Bond " << endl;
//
//			}
//            if(Schritt==Bug && thisBeadKoll) cout << "thisBeadKoll" << endl;
//            if(Schritt==Bug && nextBeadKoll) cout << "nextBeadKoll" << endl;

			if(thisBeadKoll)
			{
				if(nextBeadKoll) Chains->at(nextKollChain).EC(x-w,nextKollBead,nextKollChain,m,i,ev,d+1,thisBeadKoll,false,true);
				else {
//                    if(Schritt==Bug) cout << "else 1" << endl;

					if(rnd()<p) {
					Chains->at(nextKollChain).EC(x-w,nextKollBead,nextKollChain,m,i,ev,d+1,thisBeadKoll,true); //rechter Bead wird bewegt
					}
					else {
					Chains->at(nextKollChain).EC(x-w,nextKollBead-1,nextKollChain,m,i,ev,d+1,thisBeadKoll); //linker Bead wird bewegt
					}

				}
			}
			else{

				if(nextBeadKoll) Chains->at(nextKollChain).EC(x-w,nextKollBead,nextKollChain,m,i+rightchain,ev,d+1,thisBeadKoll,false,false);
				else {
//                    if(Schritt==Bug) cout << "else 2" << endl;
					if(rnd()<p) {
					Chains->at(nextKollChain).EC(x-w,nextKollBead,nextKollChain,m,i+rightchain,ev,d+1,thisBeadKoll,true);

					}
					else {
					Chains->at(nextKollChain).EC(x-w,nextKollBead-1,nextKollChain,m,i+rightchain,ev,d+1,thisBeadKoll);
					}
				}
			}
			return;
		}


	AnzahlFederEC++;
	if (w1<0 || w2<0) cout << "FEHLER: Wegstrecke <0!" << endl;
	if (w1<w2) {
//		cout << "Bead " << i << " wird um w1=" << w1 << "bewegt" << endl;
//		Beads.at(i).verschiebe(ev*w1);
//		i=i+1;
		Chains->at(m).EC(x-w1,i+1,m,m,i+1,ev,d+1);
		return;
		}
	if(w2<w1) {

//		cout << "Bead " << i << " wird um w2=" << w2 << "bewegt" << endl;
//		Beads.at(i).verschiebe(ev*w2);
//		i=i-1;
		Chains->at(m).EC(x-w2,i-1,m,m,i-1,ev,d+1);
		return;
		}
	if(w2==w1 && w1==x)
	{
//	    if(Schritt==Bug) cout << "Ende" << endl;

//		cout << "Bead " << i << " wird um x=" << x << "bewegt" << endl;
//		Beads.at(i).verschiebe(ev*x);
		x-=w1;
		//Chains->at(i).EC(x-w,i,m,m,i,ev);
		return;
	}

	else if(w2==w1)
	{
//		cout << "Bead " << i << " wird um w1=w2=" << w1 << "bewegt" << endl;
//		Beads.at(i).verschiebe(ev*w1);
		//x-=w1;
		if(rnd()<0.5) i++;
		else i--;
		Chains->at(m).EC(x-w1,i,m,m,i,ev,d+1);
		return;
	}


	//}

//	cout << "k= " << k << endl;

	if(x!=0) cout << "FEHLER!!! x ist nicht 0 sondern " << x << endl;
}


//void chain::EC(double x, int i,int m,int lastchain, int lastbead, vektor<double> ev, ofstream &ofbla,bool beadkoll, bool rightKollBead)
//{

//	Positionen(ofbla);
//	if(x==0) return;
//	//int i=INTrnd();//Bead i soll bewegt werde
//	if(Schritt==Bug) cout << "Bewegung von Kette " << m << " Bead " << i << " x=" << x << endl;
//	if(Schritt==Bug) cout << "lastchain =" << lastchain << " lastbead=" << lastbead << endl;
////	double phi=2*rndPI();
////	double theta=rndPI();
//	//vektor<double> ev=vektor<double>(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta));
//	//vektor<double> ev=vektor<double>(0,0,1);
//	//ev=ev/ev.Betrag();
//	//ev.Ausgabe();
//	//if(ev.Betrag()!=1) {cout << "Fehler: Betrag von ev ist " << ev.Betrag() << endl;}
//	//double x=l;

//	//while(x>0){
//	double w=x;
//	double q=x;
//
//	int nextKollChain;
//	int nextKollBead;
//	bool rightchain;
//	for(int j=0;j<M;j++)
//		{

//			for(int k=1;k<N;k++)
//				{
//					while((j==m && (k==i || k==i+1)) || (j==lastchain && k==lastbead)){//Verhindert Kollisionen mit sich selbst und letztem Bond
//						k++;
//						if(k==N) break;
//					}
//					if(k==N) break;
//

//					if(Schritt==Bug) cout << "\n Prüfe auf Kollision mit Kette " << j <<  " Beads " << k  << "," << k-1 << endl;
//					double q1;
//					double q2;
//					if(i!=0 && !(j==m && k==i-1)) q1=Kollision(Beads.at(i-1).ort, Chains->at(j).Beads.at(k-1).ort, Chains->at(j).Beads.at(k).ort, Beads.at(i).ort, ev,x);
//					else q1=x;
//					if(i!=N-1 && !(j==m && k==i+2)) q2=Kollision(Beads.at(i+1).ort, Chains->at(j).Beads.at(k-1).ort, Chains->at(j).Beads.at(k).ort, Beads.at(i).ort, ev,x);
//					else q2=x;
//					double q=min(q1,q2);
//					//cout << "---" << endl;
//					if(q<0) cout << "--------------FEHLER, q<0! q=" << q << "---------------" << endl;
//					if(q<w) //Kollision mit Chain j Beads k-1 & k
//						{
//							if(q2<q1) {

//								rightchain=true;
//							}
//							else rightchain=false;
//							w=q;
//							nextKollChain=j;
//							nextKollBead=k;
//						}
//
//				}
//		}
////	if(w!=x) w-=1e-3;// Bewegung nicht ganz bis zum Bond!
//	//cout << "x= " << x << endl;
//	double w1=ECMove(i,i+1,x,ev);
//	double w2=ECMove(i,i-1,x,ev);
////	if(Schritt==Bug) cout << "===w1= " << w1 << " ; " << " w2= " << w2 << " wKoll=" << w << "===" << endl;
//	if(w<w1 && w < w2){
////			if(Schritt==Bug) cout << "Kollision Chain " << m << " mit " << nextKollChain << endl;
//			AnzahlKollisionen++;
//			Beads.at(i).verschiebe(ev*w);
//			if(rndWAY()<0) Chains->at(nextKollChain).EC(x-w,nextKollBead,nextKollChain,m,i+rightchain,ev,ofbla);
//			else Chains->at(nextKollChain).EC(x-w,nextKollBead-1,nextKollChain,m,i+rightchain,ev,ofbla);
//			return;
//		}
//
//	AnzahlFederEC++;
//	if (w1<=0 || w2<=0) cout << "FEHLER: Wegstrecke <=0!" << endl;
//	if (w1<w2) {
////		cout << "Bead " << i << " wird um w1=" << w1 << "bewegt" << endl;
//		Beads.at(i).verschiebe(ev*w1);
//		i=i+1;
//		Chains.at(m).EC(x-w1,i,m,m,i,ev,ofbla);
//		return;
//		}
//	if(w2<w1) {
//
////		cout << "Bead " << i << " wird um w2=" << w2 << "bewegt" << endl;
//		Beads.at(i).verschiebe(ev*w2);
//		i=i-1;
//		Chains.at(m).EC(x-w2,i,m,m,i,ev,ofbla);
//		return;
//		}
//	if(w2==w1==l)
//	{
//
////		cout << "Bead " << i << " wird um x=" << x << "bewegt" << endl;
//		Beads.at(i).verschiebe(ev*l);
//		x-=w2;
//		//Chains.at(i).EC(x-w,i,m,m,i,ev);
//		return;
//	}
//
//	else if(w2==w1)
//	{
////		cout << "Bead " << i << " wird um w1=w2=" << w1 << "bewegt" << endl;
//		Beads.at(i).verschiebe(ev*w1);
//		//x-=w1;
//		if(rnd()<0.5) i++;
//		else i--;
//		Chains.at(m).EC(x-w1,i,m,m,i,ev,ofbla);
//		return;
//	}


//	//}
//
////	cout << "k= " << k << endl;

//	if(x!=0) cout << "FEHLER!!! x ist nicht 0 sondern " << x << endl;
//}

double chain::maxWeg(vektor<double> i, vektor<double> j, vektor<double> ev,double deltaE)
{
	vektor<double> t=i-j;
	vektor<double> v=ev;
	double E0=pow(t.Betrag()-x0,2);
	double B=pow(sqrt(2*deltaE/K+E0)+x0,2);
	double w=abs(-abs((t*v))+sqrt((t*v)*(t*v)-(t*t)+B));
//	cout << "w= " << w << endl;
	return w;

}

double chain::maxWegInnen(vektor<double> i, vektor<double> j,vektor<double> ev,double deltaE)
{
	vektor<double> t=i-j;
	vektor<double> v=ev;
	double E0=pow(t.Betrag()-x0,2);
	double B=pow(sqrt(2*deltaE/K+E0)-x0,2);
	double w=abs(-abs((t*v))+sqrt((t*v)*(t*v)-(t*t)+B));
//	cout << "wInnen= " << w << endl;
	return w;

}

double chain::ECMove(int i, int j, double x, vektor <double> ev)//Berechnet maximale Strecke die Bead i in Potential von Bead j zurücklegen kann
{
//	uniform_real_distribution<> distribution(0, 1);
//	auto rnd = bind(distribution, ref(*engineX));

	double deltaE=-1/Beta*log(rnd());

	if(i<0 || i>=N) cout << "Ungültiges Bead" << endl;
	if(j==-1 || j==N) return l; //Äußere Kettenglieder
//	if(Schritt==Bug) cout << "----------Berechne Weg von Kugel " << i << " in Pot von Kugel  " << j << "----------" << endl;
//	if(Schritt==Bug) cout << "deltaE init= " << deltaE << endl;
//	if(Schritt==Bug) Beads.at(i).ort.Ausgabe();
//	if(Schritt==Bug) Beads.at(j).ort.Ausgabe();
//	if(Schritt==Bug) ev.Ausgabe();
	double w=0;
	vektor<double> t=Beads.at(i).ort-Beads.at(j).ort;
//	if(Schritt==Bug) cout << "t=";
//	if(Schritt==Bug) t.Ausgabe();
//	if(Schritt==Bug) cout << "Abstand= " << t.Betrag() << endl;
	double evt=ev*t;
//	if(Schritt==Bug) cout << "ev*t=" << evt << endl;
	vektor<double> tver=t.sub(ev*(evt));
	if(evt>=0)	//Bewegung voneinander weg
	{
	//	cout << "Bewegung voneinander weg" << endl;
		if(t.Betrag()<x0)//Start in Kugel, Bewegung auf GG zu => umsonst
		{
	//		cout << "Start in Kugel" << endl;
			double q=-evt+sqrt(evt*evt-t*t+x0*x0); //Strecke bis Rand
			if (x>q) //Bewegung bis Rand möglich
			{
	//			cout << "Verlasse Kugel" << endl;
				w+=q;
			}
			else return x;
			double z=maxWeg(Beads.at(i).ort.add(ev*w),Beads.at(j).ort, ev,deltaE);
			return min(w+z,x); //Bewegung um w+z oder l

		}
		else //Start außerhalb von Kugel
		{
	//		cout << "Start außerhalb von Kugel" << endl;
			double z=maxWeg(Beads.at(i).ort,Beads.at(j).ort,ev,deltaE);
			//cout << "z=" << z << endl;
			return min(z,x);
		}
	}

	else //Bew aufeinander zu
	{
	//	cout << "Bewegung aufeinander zu" << endl;
		if (t.Betrag()<=x0)//Start in Kugel
		{
			double EMin=Federenergie(tver)-Federenergie(t);
			if (EMin<0) cout << "------------Fehler: Emin = " << EMin << "--------------" << endl;
	//		cout << "Start in Kugel" << endl;
			//double z=maxWegInnen(Beads.at(i).ort,Beads.at(j).ort,x,ev,deltaE);
			if (deltaE>=EMin) //evt*l=Entfernung bis MIN (?)
			{
				double z=-evt;
	//			cout << "Minimum wird erreicht. Strecke =" << z << " deltaE=";
				deltaE-=EMin;//Energieverbrauch bis Minimum abziehen
	//			cout << deltaE << endl;
				if(deltaE<=0) cout << "Fehler: DeltaE=" << deltaE << endl;
				if (z>=x) return x;
				z=-evt+sqrt(evt*evt-t*t+x0*x0);//Strecke bis Rand
				if (z>=x) return x;
				z+=maxWeg(Beads.at(i).ort.add(ev*w),Beads.at(j).ort,ev,deltaE);
				return min(z,x);
			}
			else return min(maxWegInnen(Beads.at(i).ort,Beads.at(j).ort,ev,deltaE),x);
		}
		else //Start außerhalb von Kugel
		{
//			if(Schritt==Bug) cout << "Start außerhalb von Kugel" << endl;

			if (tver.Betrag()>=x0) //Kugel wird nicht geschnitten
			{
//				if(Schritt==Bug) cout << "Kugel wird nicht geschnitten" << endl;
				double z=-evt; //Strecke bis min Entfernung
				if(z>=x) return x;
				else return min(z+maxWeg(Beads.at(i).ort.add(ev*z),Beads.at(j).ort,ev,deltaE),x);
			}
			else //Kugel wird mgl. geschnitten
			{
//				if(Schritt==Bug) cout << "Kugel wird möglicherweise geschnitten" << endl;
				double z=-evt-sqrt(evt*evt-t*t+x0*x0); //Strecke bis Rand
//				if(Schritt==Bug) cout << "z=" << z << endl;
				if (z>=x) return x;
//				if(Schritt==Bug) cout << "Eintritt in Kugel" << endl;
				w=maxWeg(Beads.at(i).ort.add(ev*z),Beads.at(j).ort,ev,deltaE);

				double EMin=Federenergie(tver);
				if (deltaE>=EMin)
				{
//				if(Schritt==Bug) cout << "Bewegung bis min. Entfernung " << endl;
				w=-evt; //Bew bis min Entferung
				deltaE-=EMin; //Energieverbrauch bis Minimum abziehen
				}
				else return min(w+z,x); //Delta E aufgebraucht
				if (w>=x) return x;
				//Bew bis mittlere Entfernung OK, jetzt wieder umsonst
				w=-evt+sqrt(evt*evt-t*t+x0*x0); //Entfernung bis anderer Rand
				if (w>=x) return x;
				w+=maxWeg(Beads.at(i).ort.add(ev*w),Beads.at(j).ort,ev,deltaE);
				return min(w,x);
			}
		}
	}

	cout << "Fehler" << endl;

}
