#ifndef CONTACTMAP_H
#define CONTACTMAP_H

class contactmap{
public:
    int j1;
    int j2;
    double dist;
//    int * AnzEnt;
    vector<chain>* Chains;
    bool * Contacts; //Speichert Kontakte  i1, i2, t
    contactmap(int j1, int j2, vector<chain>* Chains,double dist);
    contactmap(){};
    ~contactmap();
    void sucheKontakte(int t);



};



#endif //CONTACTMAP_H

