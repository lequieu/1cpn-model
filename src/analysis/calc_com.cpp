#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "trajectory_iterator.h"

#define TWID 500

struct AngleInfo{
    std::string vectname;
    int sitea;
    int siteb;
    double angle;
};

int main(int argc, char**argv){
    if (argc<3){
        std::cout<<"Usage: "<<argv[0]<<" <dump file> <com output file name>"<<std::endl;
        exit(1);
    }

    
    std::string dumpfilename;
    dumpfilename = argv[1];

    std::string filename;
     
    long long t;
    int nframes, natoms, dumpfreq;
    int dumpfreqprev;

    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    parser.load_dump(dumpfilename.c_str());
    std::vector<float> box_dim;
      
    //Get the number of atoms
    natoms = parser.get_numAtoms();
    box_dim = parser.get_boxDim();
    nframes = parser.get_numFrames();

    double com0[3];

    //setup file to write to
    std::ofstream ofile;
    std::string outfile;
    outfile = argv[2];
    ofile.open(outfile);
    std::vector<double> com(3,0);
    std::vector<double> comOld(3,0);

    bool firstframe = true;
    //Loop through the dump file using the parser
    for(size_t i=0; i<nframes; i++) {

        t = parser.get_current_timestep();


        //The actual functions from the parser
        parser.next_frame();
        //quats = parser.get_quat();
        //vects_f = parser.get_vect(quats,'f');
        //vects_u = parser.get_vect(quats,'u');
        //vects_v = parser.get_vect(quats,'v');
        
        com = parser.get_com(comOld);
        
        ofile << t << " " << com[0] << " " << com[1] << " " << com[2] << std::endl;
        //calculate COM
        for(size_t k=0; k<3; k++) {comOld[k] = com[k];}
        
        if (firstframe) firstframe = false;
    }  

    parser.close();
    ofile.close();

}
