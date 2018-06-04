#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "trajectory_iterator.h"

#define TWID 500

//This is a function to generate a restart frame from a specific configuration
int main(int argc, char**argv){
    
    if (argc<3){
        std::cerr<<"Usage: "<<argv[0]<<" <dump file> <frame>"<<std::endl;
        std::cerr<<"Warning: If frame equals 0, this assumes to use the last frame of a dumpfile"<<std::endl;
        exit(1);
    }

    //Initialize frame variables
    long long t; //Iterator for current timestep
    int nframes, natoms, dumpfreq, iframe;
    std::string dumpfilename;
    std::string filename;

    //Parse the inputs
    dumpfilename = argv[1];
    iframe = std::stoi(argv[2]); //Which frame is the restart point?

    
    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    parser.load_dump(dumpfilename.c_str());
    std::vector<float> box_dim;
       
    //Get the number of atoms, box dimensions, and number of frames
    natoms = parser.get_numAtoms();
    box_dim = parser.get_boxDim();
    nframes = parser.get_numFrames();

    //The number of frames to calculate with
    //Checks to see if number is greater than possible
    if (iframe > nframes) {
        std::cerr<<"Error: frame chosen is greater than number of frames"<<std::endl;
        std::cerr<<"Chosing last frame!"<<std::endl;
        iframe = nframes;
    }

    //Loop through the dump file using the parser
    for(size_t i=0; i<iframe; i++) {
        t = parser.get_current_timestep();
        parser.next_frame();
    }  

    //Set up the restart file
    std::ofstream ofile;
    ofile.open("in.restart");
    ofile<<"#Restart file for creating a new in.lammps file for a run with different purposes"<<std::endl;
    ofile<<"#<atom id >  <x>  <y>  <z>  <q1>  <q2>  <q3>  <q4>"<<std::endl;

    //Print out the positions and quaternions now
    for(size_t j=0; j<natoms; j++) {
        ofile<<j<< " "<< parser.coords_[j][0] << " "<< parser.coords_[j][1] << " "<< parser.coords_[j][2] << " ";
        ofile<<parser.quats_[j][0]<< " "<<parser.quats_[j][1]<< " "<<parser.quats_[j][2]<< " "<<parser.quats_[j][3]<<std::endl;
    }

    //Close everything
    parser.close();
    ofile.close();
}
