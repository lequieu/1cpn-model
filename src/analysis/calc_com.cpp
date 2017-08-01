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

        //MASSES are hardcoded!
    std::vector<double> masses(3+1);
    masses[1] = 196666.0000;
    masses[2] = 1950.000000;  
    masses[3] = 19500.00000;  

    std::string filename;
     
    long long t;
    int nframes, natoms, dumpfreq;
    int dumpfreqprev;

    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    parser.load_dump(dumpfilename.c_str());
    std::vector<int> atom_types;
    std::vector<float> box_dim;
    std::vector<std::vector<double>> atoms;
    std::vector<std::vector<double>> quats;
    std::vector<std::vector<double>> vects_f;
    std::vector<std::vector<double>> vects_v;
    std::vector<std::vector<double>> vects_u;
      
    //Get the number of atoms
    natoms = parser.get_numAtoms();
    atom_types = parser.get_type(); 
    box_dim = parser.get_boxDim();
    nframes = parser.get_numFrames();

    double com0[3];

    //setup file to write to
    std::ofstream ofile;
    std::string outfile;
    outfile = argv[2];
    ofile.open(outfile);

    bool firstframe = true;
    //Loop through the dump file using the parser
    for(size_t i=0; i<nframes; i++) {

        t = parser.get_current_timestep();


        //The actual functions from the parser
        atoms = parser.get_coord();
        //quats = parser.get_quat();
        //vects_f = parser.get_vect(quats,'f');
        //vects_u = parser.get_vect(quats,'u');
        //vects_v = parser.get_vect(quats,'v');
        
        //calculate COM
        double com[3], totalmass;
        int type;
        com[0] = 0.0;
        com[1] = 0.0;
        com[2] = 0.0;
        totalmass = 0;
        for (size_t iatom=0; iatom < natoms; iatom++){
          type = atom_types[iatom]-1;
          com[0] += atoms[iatom][0]*masses[type];
          com[1] += atoms[iatom][1]*masses[type];
          com[2] += atoms[iatom][2]*masses[type];
          totalmass += masses[type];
        }
        com[0] /= totalmass;
        com[1] /= totalmass;
        com[2] /= totalmass;
        
        ofile << t << " " << com[0] << " " << com[1] << " " << com[2] << std::endl;
        
        parser.next_frame();
        if (firstframe) firstframe = false;
    }  

    parser.close();
    ofile.close();

}
