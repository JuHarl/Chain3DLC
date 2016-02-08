#include "contactmap.h"

contactmap::contactmap(int j1, int j2, vector<chain>* Chains,double dist){
    this->j1=j1;
    this->j2=j2;
    this->Chains=Chains;
    this->dist=dist;
    cout << "Erstelle Contactmap Ketten" << j1 << "," << j2 << endl;
    Contacts = new bool [tMax/Schrittweite];
    for(int q=0;q<tMax/Schrittweite;q++) Contacts[q]=false;
//    for(int i=0;i<N;i++){
//        Contacts[i]=new bool * [N];
//        for(int j=0;j<N;j++){
//            Contacts[i][j]=new bool[tMax/Schrittweite];
//            for(int k=0;k<tMax/Schrittweite;k++)
//                Contacts[i][j][k]=false;
//        }
//
//    }
}

void contactmap::sucheKontakte(int t)
{

    for(int i1=0;i1<N;i1++){
        for(int i2=0;i2<N;i2++){
            if((Chains->at(j1).Beads.at(i1).ort-Chains->at(j2).Beads.at(i2).ort).Betragsquadrat()<dist){ Contacts[t]=true;
            break;
            }
            if(Contacts[t]) break;
        }
    }
}

contactmap::~contactmap(){
    cout << "Lösche Contactmap" << j1 << " " << j2 <<  endl;
//    delete[] Contacts;

}
