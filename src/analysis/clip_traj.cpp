#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "trajectory_iterator.h"



int main(int argc, char**argv){
    /*
        Script to downsample a 1cpn dump file.
        Given an initial dump file, will output a new dump file with the specified T0, TF and dT
    */
    if ((argc != 6) && (argc != 2)){
        std::cout<<"Usage #1 (print stats of file): "<<argv[0]<<" <dump file>"<<std::endl;
        std::cout<<"Usage #2 (clip file): "<<argv[0]<<" <dump file> <new dumpfile> <init time to save (T0)> <final time to save (Tf)> <frequency (dT)>"<<std::endl;
        exit(1);
    }
    bool exit_after_stats = false;
    if (argc==2) exit_after_stats = true;

   
    // Load dump
    long long ntimestep,t;
    int natoms;
    long long t0,tf,freq;
    std::ofstream ofile;

    std::string dumpfilename;
    std::string outfilename;
    dumpfilename = argv[1];
    if (!exit_after_stats){
      outfilename = argv[2];
      t0 = std::stoll(argv[3]);
      tf = std::stoll(argv[4]);
      freq = std::stoll(argv[5]);

      ofile.open(outfilename.c_str()); //make sure its blank
      ofile.close();
    }

    //Set up all vectors needed for the trajectory parser class
    TrajectoryIterator parser;
    parser.load_dump(dumpfilename.c_str());

    std::vector<float> box_dim;
    std::vector<std::vector<double>> atoms;
    std::vector<std::vector<double>> quats;
    std::vector<std::vector<double>> vects_f;
    std::vector<std::vector<double>> vects_v;
    std::vector<std::vector<double>> vects_u;
    natoms = parser.get_numAtoms();
    ntimestep = parser.get_numFrames();


    // get some features of the file
    bool firstframe = true;
    long long tstart, tend;
    int dt, dtprev, tprev;
    bool is_dt_const = true;
    for(size_t i=0; i<ntimestep; i++) {
        t = parser.get_current_timestep();
        if (firstframe) tstart = t;
        if (i==ntimestep-1) tend = t;

        if (!firstframe){
          dt = t - tprev;
          if (i==1) dtprev = dt;
          if (dtprev != dt) is_dt_const = false;
        }
        tprev = t;
        parser.next_frame();
        if (firstframe) firstframe = false;
    }  
    parser.reset();

    //Print old file stats
    std::cout << "Initial Dump File Properties: " << std::endl;
    std::cout << "Nframes="<<ntimestep<<std::endl;
    std::cout << "T0="<<tstart<< std::endl;
    std::cout << "Tf="<<tend<< std::endl;
    if (is_dt_const) std::cout << "dT="<<dt<< std::endl;
    else std::cout << "dT=notconstant!" << std::endl;
    if (exit_after_stats) exit(0); 


    //calc new file stats
    if (tf == -1) tf = tend;
    int nnewframes;
    nnewframes = (tf-t0)/freq;
 
    //Print new file stats
    std::cout << "New Dump File Properties: " << std::endl;
    std::cout << "Nframes="<<nnewframes<<std::endl;
    std::cout << "T0="<<t0<<std::endl;
    std::cout << "Tf="<<tf<<std::endl;
    std::cout << "dT="<<freq<< std::endl;

    //now clip the file
    firstframe = true;
    for(size_t i=0; i<ntimestep; i++) {
        //The actual functions from the parser
        t = parser.get_current_timestep();
         
        if ((t >= t0) && (t <= tf) && (t%freq == 0)){
          parser.append_current_frame_to_file(outfilename);
        }

        parser.next_frame();
        if (firstframe) firstframe = false;
    }  
   
}
