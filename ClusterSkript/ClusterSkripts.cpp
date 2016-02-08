#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

using namespace std;

void writeSkript(string Dateiname, string jobname, string time, string cores, string vmem, string M, string rho, string XMAX, string tMax, string Schrittweite, string Sweeps, bool CP)
{

    ofstream o;
    o.open(Dateiname);


    o << "#!/bin/bash \n\n";
    o <<" #---------- Job name\n";
    o <<"#PBS -N " << jobname << endl;

    o << "\n#--------- Mail adress\n#--------- Don't dare to use a wrong mail adress here. Two cakes!\n";
    o<<" #PBS -M julian.harland@tu-dortmund.de\n\n";

    o<<"#--------- mail for abort, begin, end :\n";
    o<<"#PBS -m abe\n\n";
    o<<"#--------- join STDOUT and STDERR\n#PBS -j oe\n\n";

    o<<"#--------- estimated (maximum) runtime\n#--------- [[hours:]minutes:]seconds[.milliseconds]\n";
    o<<"#PBS -l walltime="<< time << endl;
    o<<"#--------- MUST not be smaller than the real runtime!\n\n";

    o<<"#---------- 1 core on one node\n";
    o<<"#PBS -l nodes=1:ppn=" << cores << "\n\n";

    o<<"#---------- Maximum amount of virtual memory\n";
    o<<"#PBS -l  vmem=" << vmem << "gb\n";
    o<<"#PBS -l pvmem="<<vmem<<"gb\n\n";

    o<<"cd /data" << endl;
    o<<"FOLDER=harland/"<<jobname << endl;
    o<< "mkdir -p $FOLDER" << endl;
    o<< "cp $HOME/Chain3DLC/main $FOLDER" << endl;

    /** Hier aktivieren wenn Checkpoint eingelesen werden soll:*/
    if(CP) o<< "cp $HOME/Chain3DLC/Checkpoints/"+jobname+"/Checkpoint.txt $FOLDER" << endl;

    o<< "cd $FOLDER" << endl;
    o<< "time ./main " << M << " " << rho << " " << XMAX << " " << tMax << " " << Schrittweite << " " << Sweeps;
    if(CP) o<< " -rc";
    o<< " > stdout.log 2>&1" << endl;
    o<< "rm -f main" << endl;
    o<< "cd -" << endl;
    o<< "mv -f $FOLDER /workstation/l43/harland/cluster_output" << endl;
    o<< "date";


}

int main(int argc, char *argv[])
{
    string jobName=argv[1];
    string rho=argv[2];
    string XMAX=argv[3];
    string tMax=argv[4];
    string Schrittweite=argv[5];
    string Sweeps=argv[6];
    string time=argv[7];
    string cores=argv[8];
    string vmem=argv[9];

    int i=10;
    ofstream o;

    o.open("bash.sh");

    o << "#!/bin/bash \n\n";
    bool ReadCheckpoint=false;
    while(i<argc)
    {
        if(!strcmp( argv[i], "-rc" )) {
                ReadCheckpoint=true;
                i++;
                if(i==argc) break;
        }
        string M=argv[i];

        o << "qsub " << jobName+"rho"+rho+"M"+M+".sh" << endl;
        o << "rm -f " << jobName+"rho"+rho+"M"+M+".sh" << endl;

        writeSkript(jobName+"rho"+rho+"M"+M+".sh", jobName+"rho"+rho+"M"+M, time, cores,vmem,M,rho,XMAX,tMax,Schrittweite,Sweeps, ReadCheckpoint);
        i++;
    }




    return 0;


}
