#include "global.cpp"
int main(int argc, char *argv[])
{
//	Zufallszahlengenerator seeden

   	time_t t;
    time(&t);
    srand((unsigned int)t);
//    t=0;
    engine.seed(t);


    bool ReadCheckpoint=false;
    bool WriteCheckpoint=true;

    M=stoi(argv[1]);
    rho=stod(argv[2]);
    XMAX=stoi(argv[3]);
    tMax=stof(argv[4]);
    Schrittweite=stof(argv[5]);
    Sweeps=stof(argv[6]);
    YMAX=XMAX;
    ZMAX=XMAX;
    N=rho*XMAX*YMAX*ZMAX/M+.5;
    AnzZellen=.8*XMAX;//1.0*XMAX/(1.25*x0);
    Zellgroesse=1.0*XMAX/AnzZellen;

    string cp="Checkpoint.txt";

    for(int i=7;i<argc;i++){
        if(!strcmp( argv[i], "-rc" )) {
                ReadCheckpoint=true;
                i++;
                if(i==argc) break;
                if(strcmp(argv[i], "-nwc")) cp=argv[i];
        }
        if(!strcmp( argv[i], "-nwc" )) WriteCheckpoint=false;

    }





//	MessungEndtoEndThreads(4,pow(2,8),2,true);

	//mpq_t test;
	//mpq_init(test);
	//mpz_init(&test);
	//mpq_set_str(&test,"0.55435",10);

//	Random=new MTRand((unsigned long) time(NULL));
//	(*Random)();
//	delete Random;
    AnzahlKollisionen=0;
    AnzahlFederEC=0;
    cout << "N=" << N << " M=" << M << " V=" << XMAX << " rho=" << 1.*M*N/(XMAX*YMAX*ZMAX) << endl;
    cout << "anzZellen=" << AnzZellen << endl;
    cout << "zellgroesse=" << Zellgroesse << endl;



    list<pair<int,int>> ***Zellen;

      // Allocate memory
      Zellen = new list<pair<int,int>>**[AnzZellen];
      for (int i = 0; i < AnzZellen; ++i) {
        Zellen[i] = new list<pair<int,int>>*[AnzZellen];

        for (int j = 0; j < AnzZellen; ++j)
          Zellen[i][j] = new list<pair<int,int>>[AnzZellen];
      }
    vector<chain> * Chains;

    if(ReadCheckpoint) {
            cout << "Lese Checkpoint "<< cp << " ein" << endl;
            Chains=readCheckpoint(&engine,Zellen,cp);

    }
    else{
        Chains = new vector<chain>();

        for(int i=0;i<M;i++){
    //		cout << "Erzeuge Kette " << i << endl;
            Chains->push_back(chain(N,Chains,&engine,Zellen));

        }
	}

    if(WriteCheckpoint) cout << "Checkpointing aktiviert" << endl;

    ofstream oKonfig;
	oKonfig.open("Konfig.txt");

	oKonfig << "XMAX=" << XMAX << "\nYMAX=" << YMAX << "\nZMAX=" << ZMAX << endl;
	oKonfig << "M=" << M << "\nN=" << N << "\nl=" << l << endl;
	oKonfig << "K=" << K << "\nbeta=" << Beta << "\nx0=" << x0 << endl;
	oKonfig << "tMax=" << tMax << "\nSchrittweite= " << Schrittweite << endl;

	ofstream oPlotKonfig;
	oPlotKonfig.open("PlotKonfig.txt");
	oPlotKonfig << N << "\n" << M << "\n" <<  XMAX << "\n" << rho << "\n" << l << "\n" << K << "\n" << Beta << "\n" << x0 << "\n" << tMax << "\n" << Schrittweite<< endl;




	ofstream oGE;
	oGE.open("Gesamtenergie.txt");

	ofstream ofbla;
	ofbla.open("Pos.txt");

	ofstream ofKorr;
//	ofKorr.open("Korrelationen.txt");

	ofstream ofOrient;
	ofOrient.open("Orientierung.txt");

	ofstream ofEndToEnd;
	ofEndToEnd.open("EndToEnd.txt");

	ofstream ofHeat;
	ofHeat.open("heat01.txt");

	ofstream ofMidBeads;
	if(ReadCheckpoint){ofMidBeads.open("MidBeads.txt", ios_base::out| ios_base::app);}
	else ofMidBeads.open("MidBeads.txt");

	ofstream ofNachbar;
//	ofNachbar.open("Nachbarnt0.txt");

	ofEndToEnd << "#Schritt \t <R^2> \t Delta R^2 \n" ;

//
//	MesseOrientierung(ofOrient,Chains,ofbla,true);
//    for(int j=0;j<M;j++) {
//     Chains->at(j).NachbarDichte(1,j,ofNachbar);
//    }

//	AnzahlKollisionen=0;
//	AnzahlFederEC=0;
//    contactmap **Contactmaps = new contactmap*[M*(M-1)];
//    int k=0;
//    for(int i=0;i<M;i++){
//            for(int j=i+1;j<M;j++){
//                Contactmaps[k]=new contactmap(i,j,Chains,1);
//                k++;
//        }
//    }

//    for(int i=0;i<1e5;i++)
//        {
//        for(int k=0;k<M*N;k++){
////                cout << k << endl;
//
////            if(k%refreshrate==0) for(int j=0;j<M;j++) Chains->at(j).refreshNachbarn();
////            cout << "Resampling" << endl;
//			double phi=2*rndPI();
//			double theta=acos(1-2*rnd());
//			vektor<double> ev=vektor<double>(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta));
//	//		vektor<double> ev=vektor<double>(sin(phi),cos(phi),0);
//	//		vektor<double> ev=vektor<double>(0,1,0);
//			ev=ev/ev.Betrag();
//			int m=M*rnd();
//			int i=N*rnd();
////			if(k==0) cout << phi << theta << m << i << endl;
//			Chains->at(m).EC(l,i,m,m,i,ev);
//
//			}
//
////            cout << "Sweep beendet" << endl;
//
//        }

    int Start=Schritt;

	for(Schritt;Schritt<tMax;Schritt+=Schrittweite){
        ofNachbar.open("Nachbarn_t" + to_string(Schritt) + ".txt");
		cout << "\n!!!!!!!!!!Schritt " << Schritt <<"!!!!!!!\n" << endl;
		Positionen(ofbla,Chains);
        for(int j=0;j<M;j++) {
         Chains->at(j).NachbarDichte(1,j,ofNachbar);
        }
//        for(int j=0;j<M*(M-1)/2;j++){
////            cout << j << endl;
//            Contactmaps[j]->sucheKontakte(Schritt/Schrittweite);
//        }
        ofNachbar.close();
        if((Schritt)%(tMax/10)==0 && Schritt!=0){
//          if(Schritt+1==Bug){
            if(WriteCheckpoint){
                    cout << "Checkpoint" << endl;
                    writeCheckpoint(Chains);
                }
//            BerechneGs(Chains);

//            cout << "Berechne Korrelationen" << endl;
            //cout << "Berechne Korrelationen (g1)...";
            //ofKorr.open("Korrelationen.txt");
        	//MidBeadKorr(ofKorr,Chains);
        	//ofKorr.close();
        	//cout << "fertig, fahre fort" << endl;
        }


//        cout << rnd() << "\n" << rnd() << "\n" << rnd() << endl;

//		Positionen(ofbla,Chains);
//		LCOrientierung2(ofOrient,Chains);The file /net/l43/harland/Program…C/EC_Dense20_5/stdout.log changed on disk.
//        cout << "Berechne MidBeads...";

//		cout << "abgeschlossen!" << endl;
//        RunningStat RSEndtoEnd;
        for(int s=0;s<Schrittweite;s++){
            //if(!(Schritt==Start && ReadCheckpoint && s==0))
              //  {
            for(int j=0;j<M;j++) {
                    ofMidBeads.open("MidBeads.txt", ios_base::out| ios_base::app);
                    Chains->at(j).PosMidBead(ofMidBeads); //Berechnet Positionen der mittleren Beads
                    ofMidBeads.close();
                    //Chains->at(j).Schwerpunkt();
                }
            //}


//            cout << "EC Schritt" << endl;
		for(int k=0;k<M*N;k++){
//                cout << k << endl;

//            if(k%refreshrate==0) for(int j=0;j<M;j++) Chains->at(j).refreshNachbarn();
//            cout << "Resampling" << endl;
			double phi=2*rndPI();
			double theta=acos(1-2*rnd());
			vektor<double> ev=vektor<double>(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta));
	//		vektor<double> ev=vektor<double>(sin(phi),cos(phi),0);
	//		vektor<double> ev=vektor<double>(0,1,0);
			ev=ev/ev.Betrag();
			int m=M*rnd();
			int i=N*rnd();
//			if(k==0) cout << phi << theta << m << i << endl;
			Chains->at(m).EC(l,i,m,m,i,ev);

			}

//            cout << "Sweep beendet" << endl;

        }
//        ofEndToEnd << Schritt << "\t" << RSEndtoEnd.Mean() << "\t" << RSEndtoEnd.StandardDeviation() << endl;
//        RSEndtoEnd.Clear();
        cout << AnzahlKollisionen << " Bondkollisionen" << endl;
        cout << AnzahlFederEC << " Federkollisionen" << endl;
        cout << AnzahlEndlosschleifen << " Endlosschleifen" << endl;

//        Chains->clear();
//        delete Chains;
//        vector<chain> * Chains=readCheckpoint(&engine,Zellen);
        if((Schritt)%(tMax/10)==0 && Schritt!=0){
//          if(Schritt+1==Bug){
            cout << "Berechen Gs" << endl;
            BerechneGs(Chains);

//            cout << "Berechne Korrelationen" << endl;
            //cout << "Berechne Korrelationen (g1)...";
            //ofKorr.open("Korrelationen.txt");
        	//MidBeadKorr(ofKorr,Chains);
        	//ofKorr.close();
        	//cout << "fertig, fahre fort" << endl;
        }
		}
		cout << "Ende der Äquilibrierung" << endl;

//        for(int j=0;j<M;j++) {
//            Chains->at(j).PosMidBead(ofMidBeads); //Berechnet Positionen der mittleren Beads
//                //Chains->at(j).Schwerpunkt();
//
//        }
		//Äquilibriert (hoffentlich)

        ofMidBeads.close();

        if(WriteCheckpoint){
            cout << "Checkpoint" << endl;
            writeCheckpoint(Chains);
        }
        Schritt--;
        //ofKorr.open("Korrelationen.txt");
        //MidBeadKorr(ofKorr,Chains);
        //ofKorr.close();
        cout << "Berechne Gs" << endl;
//        BerechneGs(Chains);
//        int Entanglements[tMax/Schrittweite/10];
//        for (int i=0;i<tMax/Schrittweite/10;i++) Entanglements[i]=0;
//        for(int i=0;i<M*(M-1)/2;i++){
//            for(int t=0;t<tMax/Schrittweite;t+=10){
//                bool ent=true;
//                for(int q=0;q<10;q++){
////                    cout << "t+q=" << t+q << endl;
//                    if(Contactmaps[i]->Contacts[t+q]==0) ent=false;
//                }
////                cout << "..." << endl;
//                if(ent==true) Entanglements[t/10]++;
//            }
//
//
//        }
//        cout << "Schreibe..." << endl;

//        ofstream oEnt;
//        oEnt.open("AnzEnt.txt");
//        for(int t=0;t<tMax/Schrittweite/10;t++)
//            oEnt << t*Schrittweite*10 << "\t" <<Entanglements[t] << endl;

        RunningStat RSEnergy;

        int Genauigkeit=20;
        double maxAbstand=.5*sqrt(XMAX*XMAX+YMAX*YMAX+ZMAX*ZMAX);
        long long rValBeads[int(maxAbstand*Genauigkeit)];
        long long rValBonds[int(maxAbstand*Genauigkeit)];
        for(int i=0;i<int(maxAbstand*Genauigkeit);i++) rValBeads[i]=0;
        for(int i=0;i<int(maxAbstand*Genauigkeit);i++) rValBonds[i]=0;


        for(int q=0; q<Sweeps;q++){
            cout << "Sweep " << q << endl;
//            Chains->clear();
//            for(int i=0;i<M;i++){
    //		cout << "Erzeuge Kette " << i << endl;
//                Chains->push_back(chain(N,Chains,&engine,Zellen));

//            }
            for(int j=0;j<M;j++) ofEndToEnd << Chains->at(j).EndToEnd() << endl;
            RSEnergy.Push(Gesamtenergie(Chains));
            for(int k=0;k<M*N;k++){

//            if(k%refreshrate==0) for(int j=0;j<M;j++) Chains->at(j).refreshNachbarn();

			double phi=2*rndPI();
			double theta=acos(1-2*rnd());
			vektor<double> ev=vektor<double>(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta));
	//		vektor<double> ev=vektor<double>(sin(phi),cos(phi),0);
	//		vektor<double> ev=vektor<double>(0,1,0);
			ev=ev/ev.Betrag();
			int m=M*rnd();
			int i=N*rnd();
//			if(k==0) cout << phi << theta << m << i << endl;
			Chains->at(m).EC(l,i,m,m,i,ev);

			}
//			cout <<"Beads" << endl;
			paarVerteilungBeads(Chains,&rValBeads[0]);
//			cout << "Bonds" << endl;
			paarVerteilungBonds(Chains,&rValBonds[0]);



		}

		cout << "Energie pro Bond=" << RSEnergy.Mean() << "+-" << RSEnergy.StandardDeviation() << endl;

        ofstream o;
        o.open("paarVerteilungBonds.txt");
        for(int i=0;i<int(XMAX/2*Genauigkeit);i++)
            o << 1.*i/Genauigkeit+1./(2*Genauigkeit) << "\t" << 2.*rValBonds[i]/((N-1)*(N-1)*(M-1)*(M-1))*XMAX*YMAX*ZMAX/(4./3*PI*(pow(1.*(i+1)/Genauigkeit,3)-pow(1.*i/Genauigkeit,3)))/Sweeps << "\n" ;
        o.close();
        o.open("paarVerteilungBeads.txt");
        for(int i=0;i<int(XMAX/2*Genauigkeit);i++)
            o << 1.*i/Genauigkeit+1./(2*Genauigkeit) << "\t" << 2.*rValBeads[i]/(N*N*(M-1)*(M-1))*XMAX*YMAX*ZMAX/(4./3*PI*(pow(1.*(i+1)/Genauigkeit,3)-pow(1.*i/Genauigkeit,3)))/Sweeps << "\n" ;
        o.close();
	  // De-Allocate memory to prevent memory leak
      for (int i = 0; i < AnzZellen; ++i) {
        for (int j = 0; j < AnzZellen; ++j){
          Zellen[i][j]->clear();
          delete [] Zellen[i][j];


        }
        delete [] Zellen[i];
      }
    delete [] Zellen;
//    for(int i=0;i<M*(M-1)/2;i++) delete Contactmaps[i];
    Chains->clear();
    delete Chains;


	cout << "Ende bei Schritt" << Schritt << endl;

	return 0;
}
